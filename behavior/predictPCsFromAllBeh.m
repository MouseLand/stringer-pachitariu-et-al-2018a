function results = predictPCsFromAllBeh(dataroot, matroot, useGPU, dex)

dall=load(fullfile(dataroot, 'dbspont.mat'));

ndims0 = [1 2 3 4 8 16 32 64 128];
npc = 128;

%%
cov_res_beh = NaN*zeros(1024,length(ndims0),length(dall.db),10);
var_beh = NaN*zeros(1024,length(ndims0),length(dall.db),10);
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
		cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
		[u,s,v] = svdecon(cov);
		u = u(:,1:1024);
		v = v(:,1:1024);
		
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
			end

			cov_res_beh(:,k,d,btype) = gather_try(sout);
			var_beh(:,k,d,btype) = gather_try(vars);
			drawnow;
			hold all;
		end  
    end
    
end

%%
save(fullfile(matroot,'expv_behavior_PC.mat'),'cov_res_beh','var_beh','ndims0');


%%
