function faceTimelags(dataroot,matroot,useGPU)

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
	
	cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
	[u,s,v] = svdecon(cov);
	u = u(:,1:1024);
	v = v(:,1:1024);

    %% divide time in half
    %Ff = randn(size(Ff));
    Lblock = 60;
    fractrain = 0.5;
    [itrain, itest] = splitInterleaved(NT, Lblock, fractrain, 1);
    tic;
    if useGPU
        Ff = gpuArray(single(Ff));
	end
	
    %% use face SVDs to predict
    x = dat.beh.face.motionSVD;
    x = x - mean(x,1);
    x = x / std(x(:,1));
      
    ndims1 = ndims0(ndims0<=size(x,1));
        
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
    
    %%
    x = dat.beh.face.motionSVD;
    x = x - mean(x,1);
    x = x / std(x(:,1)) * 10;
    x    = x';
    x0 = x;    
    for k = 1:numel(tdelay)
        td = tdelay(k);
        NT = size(Ff,2);
        tpred  = max(1, 1-td) : min(NT, NT-td);
        tneur  = max(1, 1+td) : min(NT, NT+td);
        % introduce time delay
        y   = Ff(:,tneur);
        y   = y - mean(y,2);
        y   = bin2d(y, tbin, 2);
        v       = u' * y;
        
        x = x0(:,tpred);
        x = bin2d(x, tbin, 2);
        
        NT = size(y,2);
        Lblock = 60;
        fractrain = 0.5;
        
        % can increase this to average more over different splits of the data in time
        nseed = 10; 
        
        ndims1 = ndims0(ndims0<=size(x,1) & ndims0<=size(v,1));
            
        for rseed = 1:nseed
            [indtrain, indtest] = splitInterleaved(NT, Lblock, fractrain,rseed);
            clf
            itrain = find(indtrain);
            [a, b] = CanonCor2(v(:,itrain)', x(:,itrain)', .05);
            
            % prediction of neural activity
            itest  = find(indtest);
            expv = [];
            ntest = 100;
            expvPC=[];
            for kn = 1:numel(ndims1)
                n = ndims1(kn);
                clear yp;
                vp     = a(:,1:n) * b(:,1:n)' * x(:,itest);
                
                % residuals of PCs
                vres   = v(:,itest) - vp;
                expvPC(:,kn) = 1 - nanvar(vres,1,2)./nanvar(v(:,indtest),1,2);
                
                % residuals of neurons
                yp     = u * vp;% + mean(y(:,:),2);
                yres   = y(:,itest) - yp;
                expv(kn) = 1 - nanmean(nanvar(yres,1,2))/nanmean(nanvar(y(:,indtest),1,2));
            end
            
            semilogx(ndims1,expv,'*-');
            hold all;
            title([d td max(expv)]);
            axis tight;
            drawnow;
            expv_tlag(1:length(ndims1),d,k) = expv_tlag(1:length(ndims1),d,k) + expv(:)/nseed;
        end
    end
end


%%
save(fullfile(matroot,'expv_timedelay.mat'),'expv_tlag','tdelay','ndims0');















