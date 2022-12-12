function predictNeuronsFromFacePositions(dataroot,matroot)

dall=load(fullfile(dataroot, 'dbspont.mat'));

for d = [1:length(dall.db)]
    %%
    dat = load(fullfile(dataroot,sprintf('spont_%s_%s.mat',dall.db(d).mouse_name,dall.db(d).date)));
    if isfield(dat.stat, 'redcell')
        Ff = dat.Fsp(~logical([dat.stat(:).redcell]),:);
        cellpos{d} = dat.med(~logical([dat.stat(:).redcell]),:);
    else
        Ff = dat.Fsp;
        cellpos{d} = dat.med;
    end
    cellpos{d} = cellpos{d}(sum(Ff,2)>0,:);
    Ff = Ff(sum(Ff,2)>0,:);
	
    %%
    tdelay = 1;
    y    = Ff(:,(tdelay+1):end);
    [NN,NT] = size(y);
    disp('>>>>>>>>>>>>> ');
    disp([NN NT]);
    
    % 1.2 second bins
    if dat.db.nplanes == 10
        tbin = 4;
    else
        tbin = 3;
    end
    
    y = bin2d(y, tbin, 2);
    y = y - mean(y,2);
    
    [u s v] = svdecon(y);
    ncomps  = 128;
    u       = (u(:, 1:ncomps));% * s(1:ncomps,1:ncomps));
    v       = u' * y;
    
    NT = size(y,2);
    Lblock = round(60/.8);
    fractrain = 0.5;
    [indtrain, indtest] = splitInterleaved(NT, Lblock, fractrain,1);
        
    %% loop over predictors
    clf;
    x = dat.beh.face.motionSVD;
    x = x - mean(x,1);
    x = x / std(x(:,1)) * 10;
        
    x    = x(1:end-tdelay,:);
    np   = size(x,2);
    x    = squeeze(mean(reshape(x(1:floor(size(x,1)/tbin)*tbin,:),...
        tbin, [],size(x,2)),1));
    x    = reshape(x, [], np);
    x    = x';
    
    ndims0 = [1 2 3 4 8 16 32 64 128];
	
    ndims1 = ndims0(ndims0<=size(x,1) & ndims0<=size(v,1));
        
    %% low rank regression
    [a, b] = CanonCor2(v(:,indtrain)', x(:,indtrain)', .05);
    
    % prediction of neural activity
    expv = [];
    expv_neur=[];
    ntest = 100;
    for k = 1:numel(ndims1)
        clear yp;
        n       = ndims1(k);
        vp      = a(:,1:n) * b(:,1:n)' * x(:,indtest);
        
        % residuals of neurons
        yp      = u * vp;
        if n == 16
            yp0 = yp;
        end
        yres    = y(:,indtest) - yp;
        expv(k) = 1 - nanmean(nanvar(yres,1,2))/nanmean(nanvar(y(:,indtest),1,2));
        expv_neur(:,k) = 1 - nanvar(yres,1,2)./nanvar(y(:,indtest),1,2);
    end
    
    semilogx(ndims1,expv,'*-');
    hold all;
    %ylim([0 1]);
    title(max(expv));%dat{d}.db.expt_name{1});
    %ylim([0 .1])
    axis tight;
    drawnow;
    
    %expv_behavior(1:length(ndims1),d) = expv;
    expv_neurons{d} = expv_neur;
    
    %% compute 1D embedding for all datasets
%     nC=30;
%     [~,isort] = sort(u(:,1));
%     [iclust,isort] = embed1D(zscore(y(:,indtest),1,2),nC,isort);
%     yt = y(isort,indtest);
%     ytstd = max(1e-3,std(yt,1,2));
%     yt = yt ./ ytstd;
%     yp = yp0(isort,:) ./ ytstd;
%     
%     yt = my_conv2(yt,6,1);
%     yp = my_conv2(yp,1,1);
%     ccembed(d) = corr(yt(:),yp(:));
%     disp(ccembed);
%     clf;
%     subplot(2,1,1);
%     imagesc(yt,[0 .25])
%     subplot(2,1,2);
%     imagesc(yp,[0 .25])
%     drawnow;
%     ypred{d} = yp0;
%     ytest{d} = y(:,indtest);
%     isortembed{d} = isort;
    
end
%%
save(fullfile(matroot,'expv_neurons_pos.mat'),'expv_neurons','cellpos');

%%
%save('allfacepreds.mat','ypred','ytest','isortembed');












