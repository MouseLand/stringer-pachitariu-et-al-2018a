function peerPC_CCA(dataroot,matroot,useGPU,dex)

nPC = 2.^[0:10];
ndim0 =  1:512;

dall=load(fullfile(dataroot, 'dbspont.mat'));

clf;
%hold all;
lambda = 1e-5; % regularizer

%%
expvPC=[];
expvPCCov=[];
expvPCVar=[];
expv128=[];
rng('default');
for d = [1:length(dall.db)]
    %%
    dat = load(fullfile(dataroot,sprintf('spont_%s_%s.mat',dall.db(d).mouse_name,dall.db(d).date)));
    if isfield(dat.stat, 'redcell')
        Ff = dat.Fsp(~logical([dat.stat(:).redcell]),:);
        med = dat.med(~logical([dat.stat(:).redcell]),:);
    else
        Ff = dat.Fsp;
        med = dat.med;
    end
    med = med(sum(Ff,2)>0,:);
    Ff = Ff(sum(Ff,2)>0,:);
    fprintf('\nrecording %d\n',d);
    %%
    NN = size(Ff,1);
    % divide X and Y into checkerboard and use every other square
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
    
    [utrain,~,vtrain] = svdecon(Ff(ntrain,itrain));
    [utest,~,vtest] = svdecon(Ff(ntest,itrain));
    npc = min(size(vtrain,2));
    vtrain = vtrain(:,1:npc);
    vtest = vtest(:,1:npc);
    utrain = utrain(:,1:npc);
    utest = utest(:,1:npc);
    
    a = (vtrain'*vtrain + lambda*eye(npc))\ (vtest'*vtrain);
    utestrot = utest * a;
    
	
    %% PC prediction
    vtest1 = utrain' * Ff(ntrain,itest);
    vtest2 = utestrot' * Ff(ntest,itest);
    %vtest1 = vtest1 - mean(vtest1,2);
	%vtest2 = vtest2 - mean(vtest2,2);
	
    vcov = sum(vtest1.*vtest2,2);
	vcov = vcov(1:1024);
    vvar = 0.5*(sum(vtest1.^2,2) + sum(vtest2.^2,2));
    vvar = vvar(1:1024);
	
    clf
    semilogx(vcov./vvar);
    hold all;
    semilogx(cumsum(vcov)./cumsum(vvar));
    drawnow;
    
    expvPC(1:length(vcov),d) = gather_try(vcov./vvar);
    expvPCCov(1:length(vcov),d) = gather_try(vcov);
    expvPCVar(1:length(vcov),d) = gather_try(vvar);
	
    if d==dex
        exampleV1 = gather_try(vtest1);
        exampleV2 = gather_try(vtest2);
	end
	
	%% PC rank prediction
	pcrange = 1:128;
	[a, b] = CanonCor2(vtest(:,pcrange), vtrain, 1e-8);
	vtest2 = utest' * Ff(ntest,itest);	
	
	ndims1=2.^[0:7];
	expv=[];
	for k = 1:length(ndims1)
		n=ndims1(k);
		clear yp;
		vp     = a(:,1:n) * b(:,1:n)' * vtest1;
		vp     = vp(pcrange,:);
		%yp     = utest * vp - Ff(ntest,itest); 
		vres   = vtest2(pcrange,:) - vp(:,:);
		vcov = sum(vtest2(pcrange,:) .* vp, 2);
		vvar = 0.5*(sum(vtest2(pcrange,:).^2,2) + sum(vp.^2,2));
		expv(k) = gather_try(sum(vcov)/sum(vvar));
		%expv(k) = gather_try(1 - nanmean(nanvar(vres,1,2))/nanmean(nanvar(vtest2(pcrange,:),1,2)));
	end
	clf;
	plot(expv);
	title(num2str(max(expv)))
	expv128(:,d) = expv;
	drawnow;
	
end

save(fullfile(matroot,'PCpredict.mat'),'expvPC','expvPCCov','expvPCVar','expv128','ndims1','exampleV1','exampleV2');



