function smooth1Dclusters(dataroot,matroot,useGPU)
dall=load(fullfile(dataroot, 'dbspont.mat'));
dex = 1;

nC=20;
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
    
    S = zscore(Ff, 1, 2)/size(Ff,2).^.5;
    if useGPU
        S = gpuArray(single(S));
    end
    
    [u, s,v] = svdecon(S);
    U = u(:,1);
    [~,isort]=sort(U);
    [iclust,isort] = embed1D(S,nC,isort,useGPU);
    iclust = ceil(iclust / 100);
        
    %% distances
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
    if d==dex
        results.spks = Ff(isort(1:2:NN),:);
        results.beh = zscore([dat.beh.runSpeed dat.beh.pupil.area dat.beh.whisker.motionSVD(:,1)],1,1);
        Fbin = bin2d(Ff, 4, 2);
        cex = corr(Fbin');
        results.cex = cex(isort(1:2:NN),isort(1:2:NN));
        results.pos = med;
        results.iclust = iclust;
    end

end
%%

save(fullfile(matroot,'clust1D.mat'),'-v7.3','results','dex');
