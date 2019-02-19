
function [expv_all,resid,totvar,ndims0,isort,ytest,ypred] = ...
    faceAnalysis(dat, ntbin, tlag, Lam)


nt = round(120 * 30); % size of interleaved segments

y=dat.stall(:,:);
y = single(y);
ysub = mean(y,2);
y    = y - ysub;

[NN NT] = size(y);
disp('>>>>>>>>>>>>> ');
disp([NN NT]);

y = gpuArray(y);

% take SVDs of neural activity
ysm = my_conv2(y, 3, 2);
[u s v] = svdecon(ysm);
ncomps  = 128;

u       = u(:, 1:ncomps);
v       = u' * y(:,:);

%% smooth motion energy PCs
z = dat.motSVD(:,:);

nsamps = 5;

expv_all = [];
expvtrain_all = [];
ndims0 = [1 2 4 8 16 32 64 128];
resid = zeros(NN, numel(ndims0), numel(tlag));
totvar = zeros(NN, numel(ndims0), numel(tlag));


for j = 1:numel(tlag)    
    tBeh  = dat.tspont - tlag(j);
    x = interp1(dat.tVid, z, tBeh);
    X0 = gpuArray(single(x'));
      
    % compute face to neural projections from spont blocks
    resid0=[];
    totvar0=[];
    expvtrain=[];
    for rseed = 1:nsamps
		if rseed==nsamps
			[indtrain,indtest]=splitInterleaved(size(y,2),nt,.5,rseed);
		else
			[indtrain,indtest]=splitInterleaved(size(y,2),nt,.75,rseed);
		end
        
        % low rank regression
        itrain = find(indtrain);
        
        Vbin = movmean(v(:,itrain), ntbin, 2);
        Xbin = movmean(X0(:,itrain), ntbin, 2);
        
        [a, b] = CanonCor2(Vbin', Xbin', Lam);%1e3);
        a      = gather(a);
        b      = gather(b);
        
        % prediction of activity
        expv = [];
        itest = find(indtest);
        
        Ytest = movmean(y(:,itest), ntbin, 2);
        Xtest = movmean(X0(:,itest), ntbin, 2);
                
        for k = 1:numel(ndims0)
            n = ndims0(k);
            clear yp;
            vp     = a(:,1:n) * b(:,1:n)' * Xtest;
            yp     = u * vp;
            
            % residuals of neurons
            yres   = Ytest - yp;
            
            % variance explained
            expv(k) = gather(1 - nanmean(nanvar(yres,1,2))/nanmean(nanvar(Ytest,1,2)));
            resid0(:,k) = gather(nanvar(yres,1,2));
            totvar0(:,k) = gather(nanvar(Ytest,1,2));
			
			if k==5
				if tlag(j)==0.
					if rseed==nsamps
						ypred = yp;
					end
				end
			end
		end
		if rseed<nsamps
			expv_all(:,j,rseed) = expv;
			resid(:,:,j) = resid(:,:,j) + gather(resid0);
			totvar(:,:,j) = totvar(:,:,j) + gather(totvar0);
		end
	end
    
	if tlag(j) == 0.
		%%
		nC =30;
		ops.nCall = [nC 100]; % number of clusters
		ops.iPC = 1:200; % PCs to use
		ops.useGPU = 1; % whether to use GPU
		ops.upsamp = 100; % upsampling factor for the embedding position
		ops.sigUp = 1; % stddev for upsampling
		S = bin2d(y, ntbin, 2);
		S = S - my_conv2(S,100,2);
		S = zscore(S,1,2) ./ sqrt(size(S,2));

		[isort, ~, Sm] = mapTmap(S, ops);

		isort = gather_try(isort);

		ytest = gather_try(bin2d(Ytest, ntbin, 2));
		ypred = gather_try(bin2d(ypred, ntbin, 2));
	
		%%
		clf;
		subplot(2,1,1),
		imagesc(my_conv2(zscore(ytest(gather(isort), :), 1, 2),1,1),[0 .3]);
		subplot(2,1,2),
		imagesc(zscore(ypred(gather(isort), :),1,2),[0 .3]);
		cm=colormap('gray');
		colormap(flipud(cm));
		drawnow;
		pause(4);
	end
	%%
    clf;
    semilogx(ndims0,mean(expv_all(:,j,:),3),'*-');
    title([tlag(j) max(mean(expv_all(:,j,:),3))])
    
    drawnow;
    
end
resid = resid/(rseed-1);
totvar = totvar/(rseed-1);
