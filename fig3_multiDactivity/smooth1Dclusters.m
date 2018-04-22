clearvars -except dat

%% matrix for subsorting by kriging
nC = 20;
xs = 1:nC;
xn = linspace(1, nC, nC*10);

sig = 1;
d0 = (xs' - xs).^2;
d1 = (xn' - xs).^2;

K0 = exp(-d0/sig);
K1 = exp(-d1/sig);

Km = K1 / (K0 + 0.001 * eye(nC));

%%
dex = 1;
for d = 1:9
    
    Ff = dat{d}.Fsp;
     if isfield(dat{d}.stat, 'redcell')
        redcell = logical([dat{d}.stat.redcell]);
    else
        redcell = false(numel(dat{d}.stat), 1);
    end
    
    Ff = dat{d}.Fsp(~redcell(:), :);

    [NN NT] = size(Ff);
    
    %Ff = my_conv2(Ff, 1, 2);
    
    S = zscore(Ff, 1, 2)/size(Ff,2).^.5;
    S = gpuArray(single(S));
    
    [u, s,v] = svdecon(S);
    U = u(:,1);
    ncomps = 500;
    S = u(:,1:ncomps) * s(1:ncomps,1:ncomps) * v(:,1:ncomps)';
    S = gpuArray(S);
    
    %%
    [NN, NT] = size(S);
    
    iclust = zeros(NN, 1);
    [~, isort] = sort(U);
    
    nn = floor(NN/nC);
    iclust(isort) = ceil([1:NN]/nn);
    iclust(iclust>nC) = nC;
    
    sig = [linspace(3,1,25) 1*ones(1,50)];
    
    for t = 1:numel(sig)
        V = gpuArray.zeros(NT, nC, 'single');
        for j = 1:nC
            ix = iclust== j;
            V(:, j) = sum(S(ix, :),1);
        end
        
        V = my_conv2(V, sig(t), 2);
        V = normc(V);
        
        cv = S * V;
        [cmax, iclust] = max(cv, [], 2);
        
        %disp(mean(100 * cmax.^2));
        
    end
    
    
    %% distances
    med    = cell2mat({dat{d}.stat(~redcell(:)).med});
    med    = reshape(med, 2, [])';
    med    = cat(2, med*2, [dat{d}.stat(~redcell(:)).iplane]' * 35);
    
    NN = size(med,1);
    dists = zeros(NN,NN,'single');
    for j = 1:size(med,2)
        dists    = dists + (med(:,j)' - med(:,j)).^2;
    end
    dists = sqrt(dists);
    dists = dists - diag(NaN*diag(dists));
   
    for ic = 1:nC
        dc = dists(iclust==ic, iclust==ic);
        results.indist(ic,d) = nanmean(dc(:));
        results.indiststd(ic,d) = nanstd(dc(:))/sqrt(numel(dc(:))-1);
        
        dc = dists(iclust==ic, iclust~=ic);
        results.outdist(ic,d) = nanmean(dc(:));
        results.outdiststd(ic,d) = nanstd(dc(:))/sqrt(numel(dc(:))-1);
        sum(iclust==ic);
        
        results.depth{ic,d} = med(iclust==ic,3);
        results.indepth(ic,d) = nanmean(nanmean(abs(med(iclust==ic,3) - med(iclust==ic,3)')));
        results.outdepth(ic,d) = nanmean(nanmean(abs(med(iclust==ic,3) - med(iclust~=ic,3)')));
    end
    
    disp([mean(results.indist(:,d)) mean(results.outdist(:,d))])
    
    %% sorting
    [cmaxup, iclustup] = max(cv * Km', [], 2);
    100 * mean(cmaxup.^2)
    iclustup = gather(iclustup);
    [~, isort] = sort(iclustup);
        
    results.spks{d} = Ff(isort(1:2:NN),:);
    results.runspeed{d} = dat{d}.beh.runSpeed;
    if d==dex
        Fbin = bin2d(Ff, 4, 2);
        cex = corr(Fbin');
        results.cex = cex(isort(1:2:NN),isort(1:2:NN));
        results.pos = med;
        results.iclust = iclust;
    end

    %%
    clf
    plot(results.indist(:,d),results.outdist(:,d),'o')
    hold all
    plot([0 500],[0 500])
    axis tight;
    title(d);
    drawnow;
    pause(1);
    %%
    iplane = [dat{d}.stat(~redcell(:)).iplane];
    
    clf
    for ic = 1:nC
        my_subplot(6,5,ic);
        histogram(iplane(iclust==ic))
    end
    drawnow;
end
%%

save('clust1D.mat','-v7.3','results','dex');

%%

Sm = zscore(Ff,1,2) / size(Ff,2)^.5;

Sm = gpuArray(Sm(isort, :));

Sm = my_conv2(Sm, [20], [1]);
%%
clf;
imagesc(Sm, [-.2 .7]/150)
cm=colormap('gray');
colormap(flipud(cm));
%%

