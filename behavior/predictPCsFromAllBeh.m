function results = predictPCsFromAllBeh(dataroot, matroot, useGPU, dex)

dall=load(fullfile(dataroot, 'dbspont.mat'));

ndims0 = [1 2 3 4 8 16 32 64 128];
npc = 128;

clear results
%%
results.cov_res_beh = NaN*zeros(1024,length(ndims0),length(dall.db),10);
results.var_beh = NaN*zeros(1024,length(ndims0),length(dall.db),10);
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
    %% divide time in half
    %Ff = randn(size(Ff));
    Lblock = 60;
    fractrain = 0.5;
    [itrain, itest] = splitInterleaved(NT, Lblock, fractrain, 1);
    tic;
    if useGPU
        Ff = gpuArray(single(Ff));
	end
	
	cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
	[u,s,v] = svdecon(cov);
	u = u(:,1:1024);
	v = v(:,1:1024);
    %% loop over behavioral predictors
    
    for btype = [1:10]
        wmot = dat.beh.whisker.motionSVD(:,1);
        wmot = wmot * sign(mean(mean(dat.beh.whisker.motionMask(:,:,1))));
        switch btype
            %wmot = wmot * sign
            case 1
                x = zscore([dat.beh.runSpeed(:)],1,1);
            case 2
                x = zscore([dat.beh.pupil.area(:)],1,1);
            case 3
                x = zscore([dat.beh.whisker.motionSVD(:,1)],1,1);
            case 4
                x = zscore([dat.beh.runSpeed(:) dat.beh.pupil.area(:)],1,1);
            case 5
                x = zscore([dat.beh.runSpeed(:) wmot],1,1);
            case 6
                x = zscore([dat.beh.pupil.area(:) wmot],1,1);
            case 7
                x = zscore([dat.beh.runSpeed(:) dat.beh.pupil.area(:) wmot],1,1);
            case 8
                x = zscore(dat.beh.face.motionSVD(:,1),1,1);
            case 9
                x = dat.beh.face.motionSVD;
                x = x - mean(x,1);
                x = x / std(x(:,1));
            case 10
                x = dat.beh.face.motionSVD;
                x = x - mean(x,1);
                x = x / std(x(:,1));
                x = [x zscore([dat.beh.runSpeed(:) dat.beh.pupil.area(:) wmot],1,1)];
        end
      
        x    = x * 10;
        x    = x(1:end-(tdelay),:); % apply time delay
        x    = bin2d(x, tbin, 1);
        x    = x';
        
        ndims1 = ndims0(ndims0<=size(x,1));
        
		%%
		
		[atrain, btrain] = CanonCor2(Ff(ntrain,itrain)'*u, x(:,itrain)', .15);
        [atest, btest] = CanonCor2(Ff(ntest,itrain)'*v, x(:,itrain)', .15);
		k=0;
		for n = ndims1
			k=k+1;
			vp_train     = atrain(:,1:n) * btrain(:,1:n)' * x(:,itest);
			vp_test      = atest(:,1:n) * btest(:,1:n)' * x(:,itest);
        
			
			s1 = u' * Ff(ntrain,itest) - vp_train;
			s2 = v' * Ff(ntest,itest) - vp_test;
			sout = sum(s1 .* s2, 2);
			vars = sum((u' * Ff(ntrain,itest)).^2 + (v' * Ff(ntest,itest)).^2,2)/2;
			if k==6
				semilogx(sout./vars)
				ypred = [u*vp_train; v*vp_test];
				vpred = vp_test;
			end

			results.cov_res_beh(:,k,d,btype) = gather_try(sout);
			results.var_beh(:,k,d,btype) = gather_try(vars);
			drawnow;
			hold all;
		end
		
		%% sort data by 1D embedding
        if btype == 9
			nC =30;
			ops.nCall = [nC 100]; % number of clusters
			ops.iPC = 1:200; % PCs to use
			ops.useGPU = useGPU; % whether to use GPU
			ops.upsamp = 100; % upsampling factor for the embedding position
			ops.sigUp = 1; % stddev for upsampling

			[isort, ~, Sm] = mapTmap(Ff([ntrain;ntest],:),ops);
			%%
			ytest = Ff([ntrain;ntest],itest);
			ytstd = max(1e-3,std(ytest,1,2));
			yt = ytest ./ ytstd;
			yp = ypred(isort,:) ./ ytstd;

			yt = my_conv2(yt(isort,:),3,1);
			yp = my_conv2(yp,1,1);
			clf;
			imagesc(yt);

			%ccembed(d) = corr(yt(:),yp(:));
			%disp(ccembed);

			results.vtest{d} = gather_try(v'*Ff(ntest,itest));
			results.vpred{d} = gather_try(vpred);
			results.ypred{d} = gather_try(ypred);
			results.ytest{d} = gather_try(ytest);
			results.isortembed{d} = isort;
			results.xtest{d} = x(1:500,itest);
        end
		
    end
    
end

%%
results.ndims0 = ndims0;

%%
save(fullfile(matroot,'expv_behavior_PC.mat'),'-struct', 'results');

%%
