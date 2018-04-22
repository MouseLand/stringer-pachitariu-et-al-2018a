clearvars -except dat;

%%
ndims0 = [1 2 3 4 8 16 32 64 128];


tdelay = [-8*3:3:-5 -4:4 5:3:8*3];
load('../dbspont.mat');
expv_tlag = zeros(length(ndims0),length(db),length(tdelay));


clf;
%%
rng('default');
for d = [1:length(db)]
    %%
    dat = load(sprintf('../spont_%s_%s.mat',db(d).mouse_name,db(d).date));
    
    if isfield(dat.stat, 'redcell')
        Ff = dat.Fsp(~logical([dat.stat(:).redcell]),:);
    else
        Ff = dat.Fsp;
    end
    y = Ff;
    
    % 1.2 second bins
    if dat.db.nplanes == 10
        tbin = 4;
    else
        tbin = 3;
    end
    
    % bin activity to compute svd
    y = y - mean(y,2);
    y = bin2d(y, tbin, 2);
    
    [u s v] = svdecon(gpuArray(single(y)));
    ncomps  = 128;
    u       = gather(u(:, 1:ncomps));% * s(1:ncomps,1:ncomps));
    
    
    [NN NT] = size(y);
    disp('>>>>>>>>>>>>> ');
    disp([NN NT]);
    
    
    %%
    for k = 1:numel(tdelay)
        td = tdelay(k)
        NT = size(Ff,2);
        tpred  = max(1, 1-td) : min(NT, NT-td);
        tneur  = max(1, 1+td) : min(NT, NT+td);
        % introduce time delay
        y   = Ff(:,tneur);
        y   = y - mean(y,2);
        y   = bin2d(y, tbin, 2);
        v       = u' * y;
        
        x = dat.beh.face.motionSVD;
        x = x - mean(x,1);
        x = x / std(x(:,1)) * 10;
        x    = x';
        x = x(:,tpred);
        x = bin2d(x, tbin, 2);
        
        NT = size(y,2);
        Lblock = 60;
        fractrain = 0.5;
        
        nseed = 10;
        
        ndims1 = ndims0(ndims0<=size(x,1) & ndims0<=size(v,1));
            
        for rseed = 1:nseed
            [indtrain, indtest] = splitInterleaved(NT, Lblock, fractrain,rseed);
            clf
            itrain = find(indtrain);
            [a, b] = CanonCor2(gpuArray(v(:,itrain)'), gpuArray(x(:,itrain)'), .05);
            a      = gather(a);
            b      = gather(b);
            
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
            %ylim([0 1]);
            title(max(expv));%dat{d}.db.expt_name{1});
            %ylim([0 .1])
            axis tight;
            drawnow;
            expv_tlag(1:length(ndims1),d,k) = expv_tlag(1:length(ndims1),d,k) + expv(:)/nseed;
            %expvPC_tlag(:,1:length(ndims1),d,k) = expvPC;
        end
        
    end
end


%%
save('expv_timedelay.mat','expv_tlag','tdelay','ndims0');















