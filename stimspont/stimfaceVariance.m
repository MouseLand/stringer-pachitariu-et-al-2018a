function results = stimfaceVariance(dataroot,matroot,useGPU)

ndims0 = [1 2 3 4 8 16 32];

load(fullfile(dataroot,'dbstimspont.mat'));

%%
for d = 1:length(dbs)
    db = dbs(d);
    dat = load(fullfile(dataroot,...
        sprintf('stimspont_%s_%s.mat',db.mouse_name,db.date)));
    
    fprintf('recording %d\n',d);
    %%
    if isfield(dat.stat, 'redcell')
        redcell = logical([dat.stat.redcell]);
    else
        redcell = false(numel(dat.stat), 1);
    end
    gcell = ~redcell(:);
    
    %% compute stim responses
    ystim = compute_means(dat.stim.istim, dat.stim.resp,2,0);
    ystim = ystim(1:end-1,gcell,:);
    A = ystim;
    %A = A ./ max(1e-6, std(A,1,1));
    ystim1 = ystim(:,:,1) - mean(ystim(:,:,1),1);
    ystim2 = ystim(:,:,2) - mean(ystim(:,:,2),1);
    % average over neurons
    sv0= corr(ystim(:,:,1),ystim(:,:,2));
    nstim = size(ystim,1);
    results.svall{d} = 2 * mean(sum(ystim1.*ystim2,1)) ...
                        / (mean(sum(ystim1.^2,1)) + mean(sum(ystim2.^2,1)));
    disp(results.svall{d});
    results.sv{d} = diag(sv0);
    results.svraw{d} = squeeze(sum(ystim1.*ystim2,1));
    
    if d==1
        Fsp = dat.Fsp(gcell(:),:);
        [~,ix] = sort(results.sv{d},'descend');
        Fex = single(Fsp(ix(1:2000), :));
        save(fullfile(matroot,'exResponses.mat'),'Fex');
    end
    
    % stimulus variance on one repeat
    ytrain = [];
    ytest=[];
    for isti = 1:32
        isa = find(dat.stim.istim==isti);
        iss = randperm(numel(isa));
        iss = isa(iss);
        ni = numel(iss);
        ytrain = cat(1,ytrain,dat.stim.resp(iss(1:floor(ni/2)),gcell(:)));
        ytest = cat(1,ytest,dat.stim.resp(iss(floor(ni/2)+[1:floor(ni/2)]),gcell(:)));
    end
    %A = cat(3,ytrain,ytest);
    results.sv1rep{d} = diag(corr(ytrain,ytest));
    ytrain = ytrain-mean(ytrain,1);
    ytest  = ytest-mean(ytest,1);
    results.sv1repall{d} = 2 * mean(sum(ytrain .* ytest,1)) ...
                        / (mean(sum(ytrain.^2,1)) + mean(sum(ytest.^2,1)));
    disp(results.sv1repall{d})
    vnoise = var(A(:,:,1) - A(:,:,2), 1, 1) / 2;
    v1     = var(A(:,:,1), 1, 1);
    v2     = var(A(:,:,2), 1, 1);
    
    % v1 = u1 + vnoise
    % signal = (u1 + u2)/2
    % => signal = (v1 + v2 - 2*vnoise)/2
    
    results.snr{d} = (v1 + v2 - 2*vnoise) ./ (2*vnoise);

    %% timepts without behavior
    tnoface = isnan(dat.beh.face.motionSVD(:,1));
        
    % fit to periods without stim
    ttrain = ~dat.stimtpt & ~tnoface;
    
    %
    x = dat.beh.face.motionSVD(ttrain,:);
    y = dat.Fsp(gcell, ttrain);
    ysub = mean(y,2);
    y  = y - ysub;
    
    [NN NT] = size(y);
    
    % use the same binning as stimuli (3 bins)
    tbin = 3;
    x    = bin2d(x,tbin,1);
    x    = x';
    y    = bin2d(y,tbin,2);
    % time delay
    x    = x(:, 2:end);
    y    = y(:, 1:end-1);
    % normalize x
    x = x - mean(x,1);
    x = x / std(x(:,1));
    NT   = size(y,2);
    
    % divide into train and test 
    Lblock = round(75);
    fractrain = .9;
    [indtrain, indtest] = splitInterleaved(NT, Lblock, fractrain, 1);
    
    % take SVDs of neural activity
    if useGPU
        [u s v] = svdecon(gpuArray(single(y)));
    else
        [u s v] = svdecon(single(y));
    end
    ncomps  = 128;
    u       = gather_try(u(:, 1:ncomps));
    v       = u' * y;
  
    %% compute face to neural projections from spont blocks
    % low rank regression
    [a, b] = CanonCor2((v(:,indtrain)'), (x(:,indtrain))',0.5);%1e3);
    a      = gather(a);
    b      = gather(b);
    
    % prediction of activity
    expv_neurons = [];
    expv_neurons_shuffled = [];
    for n = ndims0
        clear yp;
        vp     = a(:,1:n) * b(:,1:n)' * x(:,indtest);
        yp     = u * vp;
        
        % residuals of neurons
        yres   = y(:,indtest) - yp;
        expv(n) = 1 - nanmean(nanvar(yres,1,2)) ./ nanmean(nanvar(y(:,indtest),1,2));
        expv_neurons(:,n) = 1 - (nanvar(yres,1,2)) ./ (nanvar(y(:,indtest),1,2));
        %expv_neurons_shuffled(:,n) = nanvar(y(:,indtest),1,2) - nanvar(yres,1,2);
    end
    expv = expv(ndims0);
    expv_neurons = expv_neurons(:,ndims0);
    %expv_neurons_shuffled = expv_neurons_shuffled(:,ndims0);
    clf;
    semilogx(ndims0,expv,'*-');
    
    title(max(expv));
    
    drawnow;
    
    results.facevar{d} = expv_neurons(:,ndims0==16);
    results.ev(d) = expv(ndims0==16);
   
    %% compute 1D embeddings of face and stimulus spaces
    nC = 32;
    [~,initsort] = sort(u(:,1));
    S   = zscore(u * a(:,1:10),1,2) / sqrt(32);
    [iclustup, isort] = embed1D(S, nC, initsort,useGPU);
    subplot(1,2,1),
    imagesc(S(isort,:));
    results.iface{d} = isort;
    results.sface{d} = S;
    
    S   = zscore(mean(ystim,3)',1,2) / sqrt(32);
    [iclustup, isort] = embed1D(S, nC, initsort,useGPU);
    subplot(1,2,2),
    imagesc(S(isort,:));
    results.istim{d} = isort;
    results.sstim{d} = S;
    drawnow;
    
end

%%
save(fullfile(matroot,'stimfaceRepresentation.mat'),'-struct', 'results');


