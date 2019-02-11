function results = incNeurFacePred(dataroot, matroot, useGPU)

dall=load(fullfile(dataroot, 'dbspont.mat'));

nneur0 = 2.^[8:12];

ndims0 = [1 2 4 8 16 32 64 128];
lams = [1e-3 5e-3 0.01 0.05 0.1 0.15 0.3];

npc = 1024;
npc_beh = 128;

clear results
%%

results.cov_neur = NaN*zeros(npc, length(nneur0), 10, length(dall.db));
results.var_neur = NaN*zeros(npc, length(nneur0), 10, length(dall.db));
results.cov_res_beh = NaN*zeros(npc_beh, length(ndims0), length(nneur0), length(lams), 10, length(dall.db));

rng('default');
%%
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
	
	tdelay = 1; % optimal time delay
	Ff    = Ff(:,(tdelay+1):end);
	[NN NT] = size(Ff);
	fprintf('\nrecording %d\n',d);
	disp([NN NT]);
	
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
	[NN, NT] = size(Ff);
	Ff    = squeeze(mean(reshape(Ff(:,1:floor(NT/tbin)*tbin),...
		NN, tbin, []),2));
	Ff = (Ff - mean(Ff,2));
	NT = size(Ff,2);
	
	x = dat.beh.face.motionSVD;
	x = x - mean(x,1);
	x = x / std(x(:,1));
	x    = x * 10;
	x    = x(1:end-(tdelay),:); % apply time delay
	x    = bin2d(x, tbin, 1);
	x    = x';
	
	
	%% divide time in half
	%Ff = randn(size(Ff));
	Lblock = 60;
	fractrain = 0.5;
	[itrain, itest] = splitInterleaved(NT, Lblock, fractrain, 1);
	tic;
	if useGPU
		Ff = gpuArray(single(Ff));
	end
	
	nneur1 = nneur0;
	nmax = min(length(ntest),length(ntrain));
	if nmax < nneur1(end)
		nneur1(end) = nmax;
	end
	
	for n = 1:length(nneur1)
		% do multiple subsets and avg
		for it = 1:10
			ntrain1 = ntrain(randperm(length(ntrain), nneur1(n)));
			ntest1 = ntest(randperm(length(ntest), nneur1(n)));
			[sneur, varneur, u, v] = SVCA(Ff, min(npc, nneur1(n)), ntrain1, ntest1,	itrain, itest);
			results.cov_neur(1:length(sneur),n,it,d) = gather_try(sneur);
			results.var_neur(1:length(sneur),n,it,d) = gather_try(varneur);
			
			ndims1 = ndims0(ndims0<=size(x,1));
			
			u = u(:,1:npc_beh);
			v = v(:,1:npc_beh);
			%
			for l = 1:length(lams)
				[atrain, btrain] = CanonCor2(Ff(ntrain1,itrain)'*u, x(:,itrain)', lams(l));
				[atest, btest] = CanonCor2(Ff(ntest1,itrain)'*v, x(:,itrain)', lams(l));
				k=0;
				for nd = ndims1
					k=k+1;
					vp_train     = atrain(:,1:nd) * btrain(:,1:nd)' * x(:,itest);
					vp_test      = atest(:,1:nd) * btest(:,1:nd)' * x(:,itest);
					
					
					s1 = u' * Ff(ntrain1,itest) - vp_train;
					s2 = v' * Ff(ntest1,itest) - vp_test;
					sout = sum(s1 .* s2, 2);
					vars = sum((u' * Ff(ntrain1,itest)).^2 + (v' * Ff(ntest1,itest)).^2,2)/2;
					results.cov_res_beh(:,k,n,l,it,d) = gather_try(sout);
				end
				
			end
		end
		%semilogx(results.cov_neur(:,n,it,d)./results.var_neur(:,n,it,d))
		%hold all;
		%semilogx((results.cov_neur(1:128,n,it,d) - results.cov_res_beh(:,6,n,it,d))./results.var_neur(1:128,n,it,d))
		%ylim([0 1]);
		%drawnow;
	end
	%%
	clf;
	semilogx(nanmean(results.cov_neur(:,:,:,d)./results.var_neur(:,:,:,d),3))
	hold all;
	l=7;
	semilogx(nanmean((results.cov_neur(1:128,:,:,d)-squeeze(results.cov_res_beh(:,6,:,l,:,d)))./results.var_neur(1:128,:,:,d),3))
	ylim([0 1]);
	xlim([1 npc_beh]);
	drawnow;
	
	
end

%%
results.ndims0 = ndims0;
results.lams = lams;
results.nneur0 = nneur0;

%%
save(fullfile(matroot,'increasing_neurons_SVCA_facepred.mat'),'-struct', 'results');

%%
