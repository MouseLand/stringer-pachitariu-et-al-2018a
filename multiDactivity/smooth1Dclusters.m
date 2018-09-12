function smooth1Dclusters(dataroot,matroot,useGPU, dex)
dall=load(fullfile(dataroot, 'dbspont.mat'));


% options for clustering (see activityMap.m)
nC =30;
ops.nCall = [nC 100]; % number of clusters
ops.iPC = 1:200; % PCs to use
ops.useGPU = useGPU; % whether to use GPU
ops.upsamp = 100; % upsampling factor for the embedding position
ops.sigUp = 1; % stddev for upsampling

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
    [NN,NT] = size(S);
    %
    [isort, ~, Sm] = mapTmap(S,ops);
	clf;
	imagesc(Sm,[0 2]);
	colormap('parula');
	nn = floor(NN/nC);
	iclusts = repmat([1:nC],nn,1);
	iclust = zeros(NN,1);
	iclust(isort(1:nC*nn)) = iclusts(:);
	if nC*nn < NN
		iclust(isort(nC*nn+1:end)) = nC;
	end
	    
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
		if dall.db(d).nplanes==10
			Fbin = bin2d(Ff, 4, 2);
		else
			Fbin = bin2d(Ff, 3, 2);
		end
        cex = corr(Fbin');
        results.cex = cex(isort(1:2:NN),isort(1:2:NN));
        results.pos = med;
        results.iclust = iclust;
    end

end
%%

save(fullfile(matroot,'clust1D.mat'),'-v7.3','results','dex');
