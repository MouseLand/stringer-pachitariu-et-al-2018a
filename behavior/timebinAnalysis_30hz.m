function results = timebinAnalysis_30hz(dataroot, matroot, useGPU)

load(fullfile(dataroot, 'dbspont_30hz.mat'));
dall.db = dbspont_30hz;

ndims0 = [1 2 4 8 16 32 64];
lams = [1e-4 1e-3 5e-3 0.01 0.05 0.1 0.15 0.3 0.5 .8 1.0];
tbins = [1 2 4 8 16 36 64 128 256 512 1024 2048 4096];
npc = 128;

clear results

%%
results.cov_neur = NaN*zeros(512,length(tbins),length(dall.db));
results.var_neur = NaN*zeros(512,length(tbins),length(dall.db));
results.cov_res_beh = NaN*zeros(128,length(ndims0),length(tbins),length(lams),length(dall.db));
rng('default');
%%
for d = [1:length(dall.db)]
	%%
	dat = load(fullfile(dataroot,sprintf('spont_%s_%s.mat',dall.db(d).mouse_name,dall.db(d).date)));
	
	%%
	
	if isfield(dat.stat,'redcell')
		Ff = dat.Fsp(~logical([dat.stat(:).redcell]), :);
		stat = dat.stat(~logical([dat.stat(:).redcell]));
		%med    = dat.med(~logical([dat.stat(:).redcell]),:);
	else
		Ff = dat.Fsp;
		stat = dat.stat;
		%med    = dat.med;
	end
	med = reshape([stat.med], 2, length(stat))';
	%med = med(sum(Ff,2)>0,:);
	
	Ff = Ff(sum(Ff,2)>0, :);
	
	Ff = Ff(:,~logical(dat.stimtpt));
	
	[NN, NT] = size(Ff);
	
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
	for td = 1:length(tbins)
		tbin=tbins(td);
		[NN, NT] = size(Ff);
		y    = squeeze(mean(reshape(Ff(:,1:floor(NT/tbin)*tbin),...
			NN, tbin, []),2));
		y = (y - mean(y,2));
		NT = size(y,2);
		
		tdelay = round(6/tbin); % optimal time delay
		y    = y(:,(tdelay+1):end);
		%y = y(:,~dat.stimtpt);
		[NN NT] = size(y);
		fprintf('\nrecording %d\n',d);
		disp([NN NT]);
		
		% bin face PCs in same bins
		x = dat.beh.face.motionSVD;
		x = x - mean(x,1);
		x = x / std(x(:,1));
		x    = x * 10;
		x    = x(1:end-(tdelay),:); % apply time delay
		%x= x(:,1:end/2);
		x    = bin2d(x, tbin, 1);
		x    = x';
		ndims1 = ndims0(ndims0<=size(x,1));
		
		% divide time in half
		%Ff = randn(size(Ff));
		Lblock = min(1,20*30/tbin);
		fractrain = 0.5;
		[itrain, itest] = splitInterleaved(NT, Lblock, fractrain, 1);
		tic;
		if useGPU
			y = gpuArray(single(y));
		end
		
		% SVCA
		[sneur, varneur, u, v] = SVCA(y, min(length(ntrain),length(ntest)), ntrain, ntest, itrain, itest);
		u = u(:,1:npc);
		v = v(:,1:npc);
		
		%clf;
		%subplot(1,2,1);
		semilogx(sneur./varneur)
		
		results.cov_neur(1:length(sneur),td,d) = gather_try(sneur);
		results.var_neur(1:length(sneur),td,d) = gather_try(varneur);
		hold all;
		drawnow;
		
		%%
		for l = 1:length(lams)
			[atrain, btrain] = CanonCor2(y(ntrain,itrain)'*u, x(:,itrain)', lams(l));
			[atest, btest] = CanonCor2(y(ntest,itrain)'*v, x(:,itrain)', lams(l));
			k=0;
			
			for n = ndims1
				k=k+1;
				vp_train     = atrain(:,1:n) * btrain(:,1:n)' * x(:,itest);
				vp_test      = atest(:,1:n) * btest(:,1:n)' * x(:,itest);
				
				
				s1 = u' * y(ntrain,itest) - vp_train;
				s2 = v' * y(ntest,itest) - vp_test;
				sout = sum(s1 .* s2, 2);
				%sout = max(0, sout);
				vars = sum((u' * y(ntrain,itest)).^2 + (v' * y(ntest,itest)).^2,2)/2;
				if k==length(ndims1)
					subplot(1,2,1),
					semilogx((results.cov_neur(1:npc,td,d)-sout)./vars)
					ypred = [u*vp_train; v*vp_test];
					vpred = vp_test;
					hold all;
					drawnow;
				end
				
				results.cov_res_beh(:,k,td,l,d) = gather_try(sout);
				%results.var_beh(:,k,l,d) = gather_try(vars);
				
			end
		end
		subplot(1,2,2),
		bx=cumsum(results.cov_neur(1:npc,td,d) - results.cov_res_beh(:,:,td,l,d));
		semilogx(ndims0,bx(npc,:)/sum(results.var_neur(1:npc,td,d)))
		disp(max(bx(npc,:)/sum(results.var_neur(1:npc,td,d))))
		ylim([0 max(bx(npc,:)/sum(results.var_neur(1:npc,td,d)))]);
		drawnow;
	end	
	
end

%%
results.ndims0 = ndims0;
results.lams = lams;
results.tbins = tbins;
results.npc = npc;

%%
save(fullfile(matroot,'timebins_30Hz_spont.mat'),'-struct', 'results');
