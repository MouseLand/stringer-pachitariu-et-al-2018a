function results = pcAnalysis(dataroot,matroot, useGPU)

dall=load(fullfile(dataroot, 'dbspont.mat'));
dex = 1;

results.runcorr=[];
results.indist = [];
results.outdist = [];
results.pcdist = [];
results.rcorr = [];
results.dmean = [];
results.dstd = [];
results.dnum = [];

for d = [1:length(dall.db)]
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
    fprintf('recording %d\n',d);
    % pairwise distances
    NN = size(med,1);
    dists = zeros(NN,NN,'single');
    for j = 1:3
        dists    = dists + (med(:,j)' - med(:,j)).^2;
    end
    dists = sqrt(dists);
    dists = dists - diag(NaN*diag(dists));
    
    % bin spikes and running in 1.2 second bins
    if dat.db.nplanes == 10
        tbin = 4;
    else
        tbin = 3;
    end
    
    Fbin = bin2d(Ff, tbin, 2);
    runbin = bin2d(dat.beh.runSpeed, tbin);
    pupilbin = bin2d(dat.beh.pupil.area, tbin);
    whiskbin = bin2d(dat.beh.whisker.motionSVD(:,1), tbin);
    
    %% find PCs
    if useGPU
        [u s v] = svdecon(gpuArray(single(Fbin - mean(Fbin,2))));
        u = gather(u);
        s = gather(s);
    else
        [u s v] = svdecon(single(Fbin - mean(Fbin,2)));
    end
    v = Fbin' * u;
    
    %% autocorrelation functions
    vt = Ff' * u(:,1:20);
    ac = acf_mat(vt, 400);
    results.autocorr{d} = ac(:,1:20);
    
    %% behavior correlations
    results.runcorr(d) = corr(runbin, v(:,1));
    results.pupilcorr(d) = corr(pupilbin, v(:,1));
    results.whiskcorr(d) = corr(whiskbin, v(:,1));
    
    clf
    my_subplot(2,1,1);
    hold all;
    plot(zscore(runbin))
    plot(zscore(v(:,1)));
    axis tight;
    [~,ix] = sort(u(:,1));
    
    my_subplot(2,1,2);
    imagesc(my_conv2(zscore(Fbin(ix,1:2000),1,2),[10 2],[1 2]), [0 1]);
    cm=colormap('gray');
    colormap(flipud(cm));
    
    drawnow;
    
    disp(abs([results.runcorr(d) results.pupilcorr(d) results.whiskcorr(d)]))
    
    if d == dex
        results.behavior = [runbin(:) pupilbin(:) whiskbin(:)];
        results.firstpcs  = v(:,1:3);
        results.spks = Fbin(ix,:);
    end
    
            
    %% spatial distribution of 1st PC
    rcells = zeros(NN,1);
    rcells(u(:,1)>0) = 1;
    rcells(u(:,1)<0) = 2;
    % indist: distance between cells of same sign of u(:,1)
    % indist: distance between cells of different sign of u(:,1)
    for j = 1:2
        indsj       = find(rcells==j);
        nindsj      = find(rcells~=j);
        results.indist(j,d)   = nanmean(nanmean(dists(indsj,indsj)));
        results.outdist(j,d)  = nanmean(nanmean(dists(nindsj,indsj)));
    end
    
    % change sign based on running correlation
    bcells = rcells;
    if results.runcorr(d) < 0
        bcells(rcells==1) = 2;
        bcells(rcells==2) = 1;
    end
    
    results.pcdist(d) = sum(bcells==1) / (sum(bcells==1) + sum(bcells==2));
    
    if d == dex
        results.cellmed = med;
        results.cellpc  = bcells;
    end
    
    %% distribution of correlations and repeatability
    
    % mean and standard deviation of all correlations
    ccall = corr(Fbin',Fbin');
    ccall = ccall - diag(NaN*diag(ccall));
    results.mcorr(d) = nanmean(ccall(:));
    results.stdcorr(d,:) = [prctile(ccall(:), 5) prctile(ccall(:), 95)];
        
    % split recording in half in time and compute correlations on halves
    NT = size(Fbin,2);
    Lblock = 60;
    fractrain = 0.5;
    clear indt
    [indt{1}, indt{2}] = splitInterleaved(NT, Lblock, fractrain,1);
    clear cc;
    for j = 1:2
        cc{j} = corr(Fbin(:,indt{j})');
        % set diagonal to nan's
        cc{j} = cc{j} - diag(NaN*diag(cc{j}));
    end
    cinds = ~isnan(cc{1}(:)) & ~isnan(cc{2}(:));
    % correlation of correlations
    results.rcorr(d) = corr(cc{1}(cinds), cc{2}(cinds));
    
    disp(results.rcorr(d));
    
    
    % save halves and full matrix for plotting
    if d == dex
        [~,ix] = sort(u(:,1));
        results.cc1 = cc{1}(ix(1:10:end),ix(1:10:end));
        results.cc2 = cc{2}(ix(1:10:end),ix(1:10:end));
        shifts = randi(NT-2000, NN, 1) + 1000-1;
        %%
        Fshift = zeros(size(Fbin), 'single');
        for n=1:NN
            Fshift(n,:) = circshift(Fbin(n,:), shifts(n));
        end
        
        % control: correlations of matrix shuffled in time!
        cshuff = corr(Fbin', Fshift');
        cshuff = cshuff - diag(NaN*diag(cshuff));
        dinds = find(triu(ones(NN,NN), 1));
        results.cshuff = cshuff(dinds(:));
        results.ccall = ccall(dinds(:));
    end
    
    %% correlations versus distance
    maxdist = 1.5e3;
    dbins = [0:100:maxdist];
    % only take top half of matrix
    dinds = find(triu(ones(NN,NN), 1));
    dists = dists(dinds);
    clear cc;
    cc = corr(Ff');
    cc = cc(dinds);
    
    %%
    [~,~,bin] = histcounts(dists, dbins);
    for k = 1:numel(dbins)-1
        results.dmean(k,d) = nanmean(cc(bin==k));
        results.dstd(k,d) = nanstd(cc(bin==k));
        results.dnum(k,d) = sum(bin==k & ~isnan(cc));
    end
    results.dbins =dbins;
end
    
%%
        
save(fullfile(matroot,'corr1stpc.mat'),'results','dex');
    
 
    
    
    
    
    
    
    
    
    
    
    
    
        
    
