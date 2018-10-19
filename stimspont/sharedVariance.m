function results = sharedVariance(dataroot, matroot, useGPU)

load(fullfile(dataroot,'dbstimspont.mat'));

clear results;
%
rng('default');

%%
for d = 1:length(dbs)
    db = dbs(d);
    dat = load(fullfile(dataroot,...
        sprintf('stimspont_%s_%s.mat',db.mouse_name,db.date)));
    
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
    y = dat.Fsp(gcell, ~dat.stimtpt & ~tnoface);
    
    [NN NT] = size(y);
    fprintf('recording %d\n',d);

    % bin spikes and behavior in 1.2 second bins
    tbin = 3;
    
    x    = bin2d(x, tbin, 1);
    x    = x';
    x = x - mean(x,1);
    x = x / std(x(:,1));
    
    y    = bin2d(y, tbin, 2);
        
    % subtract spont mean and divide by std before binning
    ysub = mean(y,2);
    ystd = 1e-6 + std(y,1,2);
    y    = (y - ysub)./ystd;
    
    
    % take SVD of neural activity
    if useGPU
        [u s v] = svdecon(gpuArray(single(y)));
    else
        [u s v] = svdecon(single(y));
    end
    ncomps  = 128;
    u       = gather_try(u(:, 1:ncomps));
    
    NT = size(y,2);
    Lblock = round(75);
    fractrain = 0.5;
    [indtrain, indtest] = splitInterleaved(NT, Lblock, fractrain,1);
        
    %% compute face to neural vectors during spontaneous periods
    % low rank regression
    % best with a time delay of 1 frame
    [a, b] = CanonCor2(y(:,1:end-1)', x(:,2:end)', 0.5);
    % face vectors
    uFace = normc(a(:,1:32));
    
    % spont vectors
    uSpont = u(:,1:32);
    
    %% compute stimulus vectors
    istim = dat.stim.istim(dat.stim.istim<33);
    sresp = dat.stim.resp(dat.stim.istim<33, gcell);
    
    yall = dat.Fsp(gcell,:);
    
    % subtract mean and std from spont periods
    sresp = (sresp - ysub')./ystd';
    yall  = (yall  - ysub) ./ystd;
    
    % split stimulus responses into 3 blocks
    % nstims x neurons x 3 (nstims=33 where 33 is spont)
    A = compute_means(istim, sresp, 3, 1);
    
    % stimulus vectors
    sm = mean(A(1:end,:,[2 3]),3); % mean on train blocks
    
    % stimulus projection
    sproj = normc(sm')' * yall;
    
    % face projection
    tproj = uFace' * yall;
        
    %% compute spont/face projection and compare to stims
    
    usub{1} = uFace;
    usub{2} = uSpont;
	rbin = bin2d(dat.beh.runSpeed(~tnoface,1),tbin,1);
    fbin = bin2d(dat.beh.face.motionSVD(~tnoface,1),tbin,1);
    ybin = bin2d(yall(:,~tnoface),tbin,2);
    
    for ktype = 1:2
        ist = sum(sum(isnan(A(1:end,:,:)),3),2)==0;
        Astim = A(ist,:,:);
        
        % total stimulus variance
        Vtot = sum(sum(Astim(:,:,1).*Astim(:,:,2)));
        
        % shared stim-spont subspace
        C12 = Astim(:,:, 3) * usub{ktype};
        [Ua Sa Va] = svdecon(C12);
        
        Uproj1 = normc(Astim(:,:,3)' * Ua);
        % stim-spont space 
        Uproj2 = normc(usub{ktype} * Va);
        
        uproj = usub{ktype};
        % subtract off the first dimension of stim in spont
        uproj = uproj - Uproj1(:,1) * (Uproj1(:,1)'*uproj);
        
        clear p pshared pface;
        for i = 1:2
            % stim responses in shared stim-spont without 1D subspace
            p(:,:,i) = Astim(:, :, i)  * uproj;
            % stim responses in shared stim-spont (no subtraction)
            pshared(:,:,i) = Astim(:,:,i) * Uproj2;
        end
        
        Vp = sum(sum(p(:,:,1).*p(:,:,2)));
        
        % fraction of stimulus variance in shared stim-spont space
        results.Vshared(:,d,ktype) = (sum(pshared(:,:,1).*pshared(:,:,2),1)/Vtot)';
        
        % fraction of variance in stim-spont without top PC
        results.Vall(d,ktype) = Vp/Vtot;
        
        disp([results.Vall(d,ktype) sum(results.Vshared(:,d,ktype),1)])
                
        results.Sall(:,d,ktype) = diag(Sa);
        results.Ushared{d,ktype} = Uproj2(:,1);
        
        results.runcorr(d,ktype) = corr((Uproj2(:,1)' * ybin)', rbin);
        results.facecorr(d,ktype) = corr((Uproj2(:,1)' * ybin)', fbin);
	end
	
	%% subtract shared and compute stim variance
	clear Ssub;
	results.svper = zeros(4,3);
	for ktype=1:2
		ushared = results.Ushared{d,ktype};
		Ssub{ktype} = usub{ktype} - ushared * (ushared'*usub{ktype});
	end
	ushared = results.Ushared{d,1};
	ytrain = [];
    ytest=[];
    istims=[];
    yf = zeros(32,NN);
	istimout = istim;
    for isti = 1:32
        isa = find(istim==isti);
        iss = randperm(numel(isa));
        iss = isa(iss);
        ni = numel(iss);
        ytrain = cat(1,ytrain,sresp(iss(1:floor(ni/3)),:));
        ytest = cat(1,ytest,sresp(iss(floor(ni/3)+[1:floor(ni/3)]),:));
        yf(isti,:) = mean(sresp(iss(floor(2*ni/3)+[1:floor(ni/3)]),:),1);
        istims = cat(1,istims,isti*ones(floor(ni/3),1));
		istimout(iss(floor(2*ni/3)+[1:floor(ni/3)])) = NaN;
    end
    usub{3} = normc(yf');
    Ssub{3} = usub{3} - ushared * (ushared'*usub{3});
	Ssub{3} = normc(Ssub{3});
    [u,s,v] = svdecon(double(Ssub{3}));
    %Ssub{3} = normc(u);
    %[u,s,v] = svdecon(usub{3}-mean(usub{3},2));
    %usub{3} = u;
    vrot = eye(size(v));
	%vrot = v;
	
    clf;
    for k = 1:3
        ystim = zeros(size(ytest,1),size(Ssub{k},2),2);
		ystim(:,:,1) = ytrain * Ssub{k} * vrot;
        ystim(:,:,2) = ytest * Ssub{k}  * vrot;
        yspont = y' * Ssub{k} * vrot;
		vsignal = mean((ystim(:,:,1)-mean(ystim(:,:,1),1)).*(ystim(:,:,2)-mean(ystim(:,:,2),1)));
		%vsignal = mean(ystim(:,:,1).*ystim(:,:,2),1);
		vstim      =  0.5 * (var(ystim(:,:,1), 1, 1) + var(ystim(:,:,2), 1, 1));
		vspont     = var(yspont, 1, 1);
        
        subplot(1,3,k),
        hold all;
        plot(vsignal)
        plot(vstim)
        plot(vspont)
        axis tight;
        drawnow;
        results.vsigstimspont(:,:,k,d) = [vsignal' vstim' vspont'];
	end
        
	% projections onto stims (without subtraction)
    results.istims{d} = istimout(~isnan(istimout));
    for k = 1:3
        results.projstim{d}{k} =  sresp(~isnan(istimout),:) * Ssub{k};
    end
    %results.projstim{d}{3} =  sresp(~isnan(istimout),:) * (usub{3} - results.Ushared{d,1} * (results.Ushared{d,1}'*usub{3}));
    results.projstim{d}{4} = sresp(~isnan(istimout),:) * results.Ushared{d,1};
        
	%% fit affine model to subspaces
	for k=1:3
		for additive = 0:1
			if additive == 0
				[results.fitmult(additive+1,k,d),results.multgain{d}{k},results.Rfit{d}{k}] = fitAffine(results.projstim{d}{k}, results.istims{d}, additive);
			else
				results.fitmult(additive+1,k,d) = fitAffine(results.projstim{d}{k}, results.istims{d}, additive);
				if k==3
					disp(results.fitmult(:,k,d));
				end
			end
		end
	end
	
    %% example dataset
    if d == 1
        results.sstim = normc(yf');
        results.sprojF = normc(yf')' * yall;
        
        for ktype=1:2
            results.tprojF{ktype} = usub{ktype}' * yall;
            
            % subtract shared dimension
            ushared = results.Ushared{d,ktype};
            stimsub = results.sstim - ushared * (ushared'*results.sstim);
            spontsub = usub{ktype} - ushared * (ushared'*usub{ktype});
            
            % projections onto all activity
            results.sprojS{ktype} = stimsub' * yall;
            results.tprojS{ktype} = spontsub' * yall;
            results.tshared{ktype} = ushared' * yall;
                        
		end
        
        results.runspeed = dat.beh.runSpeed;
    end
    
    
end

%%

save(fullfile(matroot,'stimvar.mat'),'-struct','results');

