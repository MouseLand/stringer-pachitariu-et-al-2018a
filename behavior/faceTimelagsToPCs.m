function faceTimelagsToPCs(dataroot,matroot,useGPU)

dall=load(fullfile(dataroot, 'dbspont.mat'));

ndims0 = [1 2 3 4 8 16 32 64 128];

tdelay = [-8*3:3:-5 -4:4 5:3:8*3];

clf;
expv_tlag = zeros(length(ndims0),length(dall.db),length(tdelay));

%%
rng('default');
for d = [1:length(dall.db)]
	%%
	dat = load(fullfile(dataroot,sprintf('spont_%s_%s.mat',dall.db(d).mouse_name,dall.db(d).date)));
	
	if isfield(dat.stat,'redcell')
		Ff = dat.Fsp(~logical([dat.stat(:).redcell]), :);
		med    = dat.med(~logical([dat.stat(:).redcell]),:);
	else
		Ff = dat.Fsp;
		med    = dat.med;
	end
	med = med(sum(Ff,2)>0,:);
	Ff = Ff(sum(Ff,2)>0,:);
	
	[NN NT] = size(Ff);
	fprintf('\nrecording %d\n',d);
	disp([NN NT]);
	Ff0 = Ff;
	%% divide X and Y into checkerboard and use every other square
	y = round(med(:,1));
	ymax=max(med);
	ymax = ymax(1);
	nby = floor(ymax / 16);
	ytrain = ([1:2:16]-1) * nby + [1:nby-10]';
	ytrain = ytrain(:)';
	ytrain = repmat(ytrain,NN,1);
	nt= size(ytrain,2);
	ntrain = find(sum(repmat(y,1,nt) == ytrain, 2)>0);
	ytest = ([2:2:16]-1) * nby + [1:nby-10]';
	ytest = ytest(:)';
	ytest = repmat(ytest,NN,1);
	nt = size(ytest,2);
	ntest = find(sum(repmat(y,1,nt) == ytest, 2)>0);
	
	%% bin spikes in 1.2 second bins
	if dat.db.nplanes==10
		tbin = 4;
	else
		tbin=3;
	end
	if useGPU
		Ff0 = gpuArray(single(Ff0));
	end
	Ff0 = Ff0 - mean(Ff0,2);
	NT = size(Ff0,2);
	
	% use face SVDs to predict
	x = dat.beh.face.motionSVD;
	x = x - mean(x,1);
	x = x / std(x(:,1));
	x = x*10;
	x0 = x';
	if useGPU
		x0 = gpuArray(single(x0));
	end
	
	%% divide time in half
	clf;
	%%
	for k = 1:numel(tdelay)
		td = tdelay(k);
		tpred  = max(1, 1-td) : min(NT, NT-td);
		tneur  = max(1, 1+td) : min(NT, NT+td);

		Ff = bin2d(Ff0(:,tneur), tbin, 2);
		Ff = Ff - mean(Ff,2);
		x  = bin2d(x0(:,tpred), tbin, 2);

		Lblock = 60;
		fractrain = 0.5;
		[itrain, itest] = splitInterleaved(size(Ff,2), Lblock, fractrain, 1);

		cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
		if useGPU
			[u,s,v] = svdecon(gpuArray(single(cov)));
		else
			[u,s,v] = svdecon(cov);
		end
		u = u(:,1:1024);
		v = v(:,1:1024);


		s1 = u' * Ff(ntrain, itest);
		s2 = v' * Ff(ntest, itest);
		sneur = sum(s1 .* s2, 2);
		varneur = sum(s1.^2 + s2.^2,2)/2;

		[atrain, btrain] = CanonCor2(Ff(ntrain, itrain)'*u, x(:,itrain)', .15);
		[atest, btest] = CanonCor2(Ff(ntest, itrain)'*v, x(:,itrain)', .15);

		n=16;
		vp_train     = atrain(:,1:n) * btrain(:,1:n)' * x(:,itest);
		vp_test      = atest(:,1:n) * btest(:,1:n)' * x(:,itest);

		s1 = s1 - vp_train;
		s2 = s2 - vp_test;
		sout = sum(s1 .* s2, 2);

		if k==13
			subplot(1,2,1),
			semilogx((sneur-sout)./varneur,'k','linewidth',2)
		else
			subplot(1,2,1),
			semilogx((sneur-sout)./varneur)
			hold all;
		end

		if k>1
			tpr = tt;
		end
		tt = cumsum(sneur-sout) ./ cumsum(varneur);
		tt = tt(128,:);
		subplot(1,2,2);
		if k>1
			plot(tdelay(k-1:k),[tpr tt],'-*');
		else
			plot(tdelay(k),[tt],'*');
		end

		tlag_cov_res_beh(:,k,d) = gather_try(sout);
		tlag_var_beh(:,k,d) = gather_try(varneur);
		tlag_cov_neur(:,k,d) = gather_try(sneur);
		title(d);
		drawnow;
		hold all;
	end
end

%%
save(fullfile(matroot,'expv_timedelay.mat'),'tlag_cov_neur','tlag_cov_res_beh','tlag_var_beh','tdelay','ndims0');















