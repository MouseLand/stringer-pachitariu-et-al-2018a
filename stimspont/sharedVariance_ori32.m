function results = sharedVariance_ori32(dataroot, matroot, useGPU)

%%

load(fullfile(dataroot,'dbori32.mat'));
% first orientation dataset doesn't have spont periods
dbs = dbs(2:end);

clear results;
%
rng('default');

%%
lam = [.05 1 0.5 0.1];

%%
for d = 1:length(dbs)
	
    db = dbs(d);
    dat = load(fullfile(dataroot,...
        sprintf('orispont_%s_%s.mat',db.mouse_name,db.date)));
    
    %
    if isfield(dat.stat, 'redcell')
        redcell = logical([dat.stat.redcell]);
    else
        redcell = false(numel(dat.stat), 1);
    end
    gcell = ~redcell(:);
    
    % spontaneous activity data
    x = dat.beh.face.motionSVD;
    
    % some videos may not have captured - if so there are nan's
    tnoface = isnan(x(:,1));
    
    % spontaneous activity
    % ~dat.stimtpt
    x = x(~dat.stimtpt & ~tnoface,:);
	fpc = size(x,2);
    y = dat.Fsp(gcell, ~dat.stimtpt & ~tnoface);
    
    [NN NT] = size(y);
    fprintf('recording %d\n',d);

    % bin spikes and behavior in 1.2 second bins
    tbin = 3;
    
    x    = bin2d(x, tbin, 1);
    x    = x';
    x    = x - mean(x,1);
    x     = x / std(x(:,1));
	%x     = x*10;
    
    y    = bin2d(y, tbin, 2);
    
    % subtract spont mean and divide by std before binning
    ysub = mean(y,2);
    ystd = 1e-6 + std(y,1,2);
    y    = (y - ysub)./ystd;
    
    NT = size(y,2);
    Lblock = round(60);
    fractrain = 0.7;
    [indtrain, indtest] = splitInterleaved(NT, Lblock, fractrain,2);
        
    %% compute face to neural vectors during spontaneous periods
    % low rank regression
    % best with a time delay of 1 frame
	ytrain = y(:,indtrain);
	ytrain = ytrain(:,1:end);
	ytrain = ytrain - mean(ytrain,2);
	xtrain = x(:,indtrain);
	xtrain = xtrain(:,1:end);
	ytest = y(:,indtest);
	ytest = ytest(:,1:end);
	ytest = ytest - mean(ytest,2);
	xtest = x(:,indtest);
	xtest = xtest(:,1:end);
    
	% spont vectors
    % take SVD of neural activity
    if useGPU
        [u s v] = svdecon(gpuArray(single(ytrain)));
    else
        [u s v] = svdecon(single(ytrain));
    end
    ncomps  = 128;
	ncomps  = min(ncomps,size(u,2));
    u       = gather_try(u(:, 1:ncomps));
	v       = gather_try(v(:, 1:ncomps)*s(1:ncomps,1:ncomps));
    
	%
	uSpont = normc(u(:,1:min(ncomps,size(u,2))));
	%
	[a, b] = CanonCor2(ytrain', xtrain', lam(d));
    % face vectors
    uFace = normc(a(:,1:32));
	ypred = a(:,1:32)*b(:,1:32)'*xtest;
	fprintf('face pred: %2.3f\n',1-mean(sum((ypred-ytest).^2,2))/mean(sum(ytest.^2,2)));
    
    % compute face PCs in stimulus bins
	xst = dat.beh.face.motionSVD_stim;
    xst = xst'; 
	xst = xst - nanmean(xst,1);
	xst = xst / nanstd(xst(:,1));
	
	
	rsp = dat.beh.runSpeed_stim;
	rsp = rsp - nanmean(rsp);
	
	% compute stimulus vectors
    istim = dat.stim.istim(dat.stim.istim<33);
    sresp = dat.stim.resp(dat.stim.istim<33, gcell);
	rsp   = rsp(dat.stim.istim<33);
	xst   = xst(:, dat.stim.istim<33);
	
    
    yall = dat.Fsp(gcell,:);
    
    % subtract mean and std from spont periods
    sresp = (sresp - ysub')./ystd';
    yall  = (yall  - ysub) ./ystd;
    
        
    % compute spont/face projection and compare to stims
    
    usub{1} = uFace;
    usub{2} = uSpont;
	rbin = bin2d(dat.beh.runSpeed(~tnoface,1),tbin,1);
    fbin = bin2d(dat.beh.face.motionSVD(~tnoface,1),tbin,1);
    ybin = bin2d(yall(:,~tnoface),tbin,2);

	yresp_train = [];
    yresp_test=[];
    istims=[];
    yf = zeros(32,NN);
	ystim = zeros(32,NN,3);
	istimout = istim;
    for isti = 1:32
        isa = find(istim==isti);
        iss = randperm(numel(isa));
        iss = isa(iss);
        ni = numel(iss);
        yresp_train = cat(1,yresp_train,sresp(iss(1:floor(ni/3)),:));
        yresp_test = cat(1,yresp_test,sresp(iss(floor(ni/3)+[1:floor(ni/3)]),:));
		ystim(isti,:,1) = mean(sresp(iss(floor(ni/3)+[1:floor(ni/3)]),:),1);
		ystim(isti,:,2) = mean(sresp(iss(floor(ni/3)+[1:floor(ni/3)]),:),1);
        ystim(isti,:,3) = mean(sresp(iss(floor(2*ni/3)+[1:floor(ni/3)]),:),1);
        istims = cat(1,istims,isti*ones(floor(ni/3),1));
		istimout(iss(floor(2*ni/3)+[1:floor(ni/3)])) = NaN;
	end
	results.facest{d} = xst(:,~isnan(istimout));
	results.rsp{d} = rsp(~isnan(istimout));
	results.istims{d} = istimout(~isnan(istimout));
    
    ystim_avg = ystim(:,:,3);
	usub{3} = normc(ystim_avg');
	
	%%	
	clear Ssub;
    for k = 1:2		
        % total stimulus variance
        %Vtot = sum(sum(ystim(:,:,1).* ystim(:,:,2)));
		%Vtot = sum(sum(yresp_train.^2 + yresp_test.^2)/2);
		Vtot = sum(sum(yresp_train .* yresp_test));
        
		uproj = usub{k};
        % shared stim-spont subspace
        C12 = ystim_avg * uproj;
        [Ua Sa Va] = svdecon(C12);
        
        %Uproj1 = normc(ystim_test' * Ua);
		Uproj1 = normc(ystim_avg' * Ua);
        % stim-spont space 
        Uproj2 = normc(uproj * Va);
        
        % subtract off the first dimension of stim in spont
        uproj = normc(uproj - Uproj2(:,1) * (Uproj2(:,1)'*uproj));
        %Ssub{k} = usub{k} - results.Ushared{d,k} * (results.Ushared{d,k}'*usub{k});
		Ssub{k} = normc(uproj);
		
        clear p pshared pface;
        for i = 1:2
            % stim responses in shared stim-spont without 1D subspace
			if i == 1
				p(:,:,i) = yresp_train  * Uproj2(:,2:end);
				% stim responses in shared stim-spont (no subtraction)
				pshared(:,:,i) = yresp_train * Uproj2(:,1:end-1);
			else
				p(:,:,i) = yresp_test  * Uproj2(:,2:end);
				% stim responses in shared stim-spont (no subtraction)
				pshared(:,:,i) = yresp_test * Uproj2(:,1:end-1);
			end
		end
        
		%p = p - mean(p,1);
		%pshared = pshared - mean(pshared,1);
        Vp = sum(sum(p(:,:,1).*p(:,:,2)));
        		
        % fraction of stimulus variance in shared stim-spont space
        results.Vshared(:,d,k) = (sum(pshared(:,:,1).*pshared(:,:,2),1)/Vtot)';
        
        % fraction of variance in stim-spont without top PC
        results.Vall(d,k) = Vp/Vtot;
        
        disp([results.Vall(d,k) sum(results.Vshared(:,d,k),1)])
                
        results.Sall(:,d,k) = diag(Sa);
        results.Ushared{d,k} = Uproj1(:,1);
        
        results.runcorr(d,k) = corr((Uproj1(:,1)' * ybin)', rbin);
        results.facecorr(d,k) = corr((Uproj1(:,1)' * ybin)', fbin);
		
		% correlation between stim-spont shared and 1st spont PC
		results.shared_corr_spont_spont(d,k) = corr((Uproj1(:,1)' * ytest)', (usub{1}(:,1)' * ytest)');
		results.shared_corr_face_spont(d,k) = corr((Uproj1(:,1)' * ytest)', (usub{2}(:,1)' * ytest)');
		
	end
	
	%% compute face related variance
	results.svper = zeros(4,3);
    
	for k = 1:2
		vexp = 1 - mean(mean((Ssub{k}'*ytest - Ssub{k}'*a(:,1:32)*b(:,1:32)'*xtest).^2,2)) /...
			mean(var(Ssub{k}'*ytest, 1, 2));
		results.face1st_r2(k,d) = corr((Ssub{k}(:,1)'*ytest)', (Ssub{k}(:,1)'*a(:,1:32)*b(:,1:32)'*xtest)');
		results.facepred_vexp(k,d) = vexp;
		results.facepred_var(k,d) = sum(var(Ssub{k}'*ytest, 1, 2) - mean((Ssub{k}'*ytest - Ssub{k}'*a(:,1:32)*b(:,1:32)'*xtest).^2, 2));
		if k == 1
			disp([vexp, results.facepred_var(k,d)]);
		end
	end
	%
	
	results.sigvar_avg(d) = mean(mean((yresp_train-mean(yresp_train,1)) .* ...
		                          (yresp_test-mean(yresp_train,1)),1))/...
		                           mean(var(yresp_train,1,1));

	% face-related variance in the shared subspace
    ushared = results.Ushared{d,1};
	yt = ushared' * ytrain;
	as = (xtrain*xtrain' + 1000*eye(fpc)) \ (xtrain*yt(:));
	yp = xtest'*as;
	cc = corr(yp(:), ytest' * ushared)
	results.faceushared_r2(d) = cc;
						
	%% compute stim-only space
	Ssub{3} = usub{3} - ushared * (ushared'*usub{3});
	Ssub{3} = normc(Ssub{3});
	Ssub{4} = normc(ushared);
	%%
    clf;
    for k = 1:24	
		if k < 5
			uproj = Ssub{k};
		else
			uproj = normc(randn(NN,32));
		end
        yst = zeros(size(yresp_test,1),size(uproj,2),2);
		yst(:,:,1) = yresp_train * uproj;
        yst(:,:,2) = yresp_test * uproj;
        yspont = ytest' * uproj;
		%vsignal = mean((yst(:,:,1) - mean(yst(:,:,1),1)).*(yst(:,:,2) - mean(yst(:,:,1),1)));
		%vsignal = mean(ystim(:,:,1).*ystim(:,:,2),1);
		vstim      =  0.5 * (var(yst(:,:,1), 1, 1) + var(yst(:,:,2), 1, 1));
		vspont     = var(yspont, 1, 1);
		vnoise     = 0.5*var(yst(:,:,1) - yst(:,:,2), 1, 1);
        vsignal    = vstim - vnoise;
		vsnr       = sum(vsignal) ./ sum(vnoise);
		
		if k < 5
			subplot(1,4,k),
			hold all;
			plot([vsignal' vstim' vspont'])
			axis tight;
			drawnow;
		end
        results.vsigstimspont(:,k,d) = [sum(vsignal) sum(vstim) sum(vspont) (vsnr)];
				
		% decode direction
		[results.decoding(k,d,1),~] =  decoder(istims, istims, ...
			yst(:,1:min(size(yst,2),32),1), yst(:,1:min(size(yst,2),32),2), 1);
		
		% decode orientation
		istims_ori = istims;
		istims_ori(istims_ori>16) = istims_ori(istims_ori>16) - 16;
		[results.decoding(k,d,2),~] =  decoder(istims_ori, istims_ori, ...
			yst(:,1:min(size(yst,2),32),1), yst(:,1:min(size(yst,2),32),2), 2);
		
		%results.decoding(k,d) =  decoder(istims, istims, v, v2);
		title(results.decoding(k,d,:));
		drawnow;
	end
        
	% projections onto stims (without subtraction)
	for k = 1:3
        results.projstim{d}{k} =  sresp(~isnan(istimout),:) * Ssub{k};
    end
    %results.projstim{d}{3} =  sresp(~isnan(istimout),:) * (usub{3} - results.Ushared{d,1} * (results.Ushared{d,1}'*usub{3}));
    results.projstim{d}{4} = sresp(~isnan(istimout),:) * results.Ushared{d,1};
	results.projstim{d}{5} = sresp(~isnan(istimout),:) * results.Ushared{d,2};
	for k = 1:2
		results.shared_corr_spont_stim(d,k) = corr(results.projstim{d}{3+k}, ...
			                                       sresp(~isnan(istimout),:) * usub{1}(:,1));
	    results.shared_corr_face_stim(d,k) = corr(results.projstim{d}{3+k}, ...
			                                       sresp(~isnan(istimout),:) * usub{2}(:,1));
	end
		
        
	%% fit affine model to subspaces
	for k=1:3
		for additive = 0:1
			if additive == 0
				[results.fitmult(additive+1,k,d),results.multgain{d}{k},results.Rfit{d}{k},sm] = fitAffine(results.projstim{d}{k}(:,1:32), results.istims{d}, additive);
			else
				results.fitmult(additive+1,k,d) = fitAffine(results.projstim{d}{k}(:,1:32), results.istims{d}, additive);
				if k==3
					disp(results.fitmult(:,k,d));
				end
			end
		end
	end
	
	%% predict the multiplicative gain with the face
	lamg=[100 100 100 100];
	%for rg=1:10
	[ttrain,ttest]=splitInterleaved(length(results.multgain{d}{3}),2,0.8,1);
	mg = results.multgain{d}{3};
	fg = results.facest{d};
	mg = mg-mean(mg);
	fg = fg-mean(fg,2);
	
	a=(fg(:,ttrain)*fg(:,ttrain)' + lamg(d)*eye(fpc))\(fg(:,ttrain)*mg(ttrain));
	
	r2train = corr(fg(:,ttrain)'*a, mg(ttrain));
	clf;
	plot(fg(:,ttest)'*a);
	hold all;
	plot(mg(ttest))
	r2test = corr(fg(:,ttest)'*a, mg(ttest))
	results.face_pred_gain_r2(d) = r2test;
	results.face_shared_gain_r2(d) = corr(results.projstim{d}{4}(ttest), mg(ttest));
	results.spont_shared_gain_r2(d) = corr(results.projstim{d}{5}(ttest), mg(ttest));
	%end
	
    %% example dataset
    if d == 4
        results.sstim = normc(ystim_avg');
        results.sprojF = normc(ystim_avg')' * yall;
        
        for k=1:2
            results.tprojF{k} = usub{k}' * yall;
            
            % subtract shared dimension
            ushared = results.Ushared{d,k};
            stimsub = Ssub{3};
            spontsub = Ssub{k};
            
            % projections onto all activity
            results.sprojS{k} = stimsub' * yall;
            results.tprojS{k} = spontsub' * yall;
            results.tshared{k} = ushared' * yall;
                        
		end
        
        results.runspeed = dat.beh.runSpeed;
    end
    
    
end

%%

save(fullfile(matroot,'stimvar_ori32.mat'),'-struct','results');

