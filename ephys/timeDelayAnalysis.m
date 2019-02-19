function [cov_neur_t, var_neur, cov_res_beh_t, acf_beh] = timeDelayAnalysis(dat,tlag0)
useGPU=1;
%%
ndims0 = [1 2 4 8 16 32 64 128];
lams = [1e-3 0.01 0.05 0.15 0.3 .8 1.0];

npc = 128;

results.cov_neur_t = NaN*zeros(npc,length(tlag0));
results.var_neur = NaN*zeros(npc,1);
%results.cov_res_beh = NaN*zeros(npc,length(ndims0),length(tbins),length(lams));
results.cov_res_beh_t = NaN*zeros(npc,length(ndims0),length(tlag0),length(lams));

rng('default');

%% divide probes and use every other segment for ntrain and ntest
NN = length(dat.Wh);
ntrain = [];
ntest = [];
ymax = max(dat.Wh);
nby = floor(ymax / 16);
%
for k =1:length(unique(dat.iprobe))
	ineu = find(dat.iprobe==k);
	ineu = ineu(:);
	NN = length(ineu);
	ytrain = ([1:4:16]-1) * nby + [1:nby-10]';
	ytrain = ytrain(:)';
	ytrain = repmat(ytrain,NN,1);
	nt= size(ytrain,2);
	ytest = ([3:4:16]-1) * nby + [1:nby-10]';
	ytest = ytest(:)';
	ytest = repmat(ytest,NN,1);
	ntrain = [ntrain; ineu(find(sum(repmat(dat.Wh(ineu),1,nt) == ytrain, 2)>0))];
	ntest = [ntest; ineu(find(sum(repmat(dat.Wh(ineu),1,nt) == ytest, 2)>0))];
end

%%

Ff = dat.stall(:,:);

[NN NT] = size(Ff);
disp('>>>>>>>>>>>>> ');
disp([NN NT]);

z = dat.motSVD(:,:);

%acf_beh = acf(z,30*8);
acf_beh = 0; 

ndims1 = [1 2 4 8 16 32 64];

%% bin spikes in 30 ms bins
tbin = 1;
td = 1;
y    = bin2d(Ff, tbin, 2);
y = (y - mean(y,2));

% divide time in half
[NN, NT] = size(y);
Lblock = max(1,60*30/tbin);
fractrain = 0.5;
[itrain, itest] = splitInterleaved(NT, Lblock, fractrain, 1);
tic;


tlag=tlag0;
itrain(1 : tlag(end)*15) = 0;
itrain(end - tlag(end)*15 : end) = 0;
itest(1 : tlag(end)*15) = 0;
itest(end - tlag(end)*15 : end) = 0;

% SVCA
[sneur, varneur, u, v] = SVCA(y, min(1024, min(length(ntrain),length(ntest))), ntrain, ntest, itrain, itest);
u = u(:,1:npc);
v = v(:,1:npc);

subplot(1,2,1);
semilogx(sneur./varneur);
results.cov_neur = gather_try(sneur(1:npc));
results.var_neur(1:npc) = gather_try(varneur(1:npc));
hold all;
drawnow;

%%
[u1,~,v1] = svdecon(y);
[~,isort]=sort(u1(:,1));

%%

for j = 1:numel(tlag)
	clf;
	
	
	tshift = tlag(j) * 15;
	
	%% shift svca
	
	if tshift < 0
		tshift=abs(tshift);
		itrain1 = find(itrain) - tshift;
		itest1 = find(itest) - tshift;
		itrain2 = find(itrain) + tshift;
		itest2 = find(itest) + tshift;
	else
		itrain1 = find(itrain) + tshift;
		itest1 = find(itest) + tshift;
		itrain2 = find(itrain) - tshift;
		itest2 = find(itest) - tshift;		
	end
	cov = y(ntrain,itrain1) * y(ntest,itrain2)';
	[u0,~,v0] = svdecon(cov);
	u0 = u0(:,1:npc);
	v0 = v0(:,1:npc);
	s1 = u0' * y(ntrain,itest1);
	s2 = v0' * y(ntest,itest2);
	sneur = sum(s1 .* s2, 2);
	varneur = sum(s1.^2 + s2.^2,2)/2;
	
	results.cov_neur_t(1:length(sneur),j)=gather_try(sneur);
	
	
	% shift behavior
	tBeh  = dat.tspont - tlag(j);
	motSVD = interp1(dat.tVid, z, tBeh);
	% bin face PCs in same bins after shifting
	x    = motSVD;
	x    = x - mean(x,1);
	x    = x / std(x(:,1));
	x    = x * 10;
	x    = bin2d(x, tbin, 1);
	x    = x';
	
	[NN NT] = size(y);	
	if useGPU
		y = gpuArray(single(y));
	end
	
	%%
	for l = 1:length(lams)
		[atrain, btrain] = CanonCor2(y(ntrain,itrain)'*u, x(:,itrain)', lams(l));
		[atest, btest] = CanonCor2(y(ntest,itrain)'*v, x(:,itrain)', lams(l));
		k=0;
		
		xtrain = btrain' * x(:,itest);
		xtest = btest' * x(:,itest);
		ytrain =  atrain' * u' * y(ntrain,itest);
		ytest =   atest' * v' * y(ntest,itest);
		
		%cc = 0.5 * (ccf(ytrain(1:16,:)', xtrain(1:16,:)', 1000) + ...
		%	ccf(ytest(1:16,:)', xtest(1:16,:)', 1000));
		
		%ccA(:,:,td,l,d) = cc;
		
		for n = ndims1
			k=k+1;
			vp_train     = atrain(:,1:n) * btrain(:,1:n)' * x(:,itest);
			vp_test      = atest(:,1:n) * btest(:,1:n)' * x(:,itest);
			
			
			s1 = u' * y(ntrain,itest) - vp_train;
			s2 = v' * y(ntest,itest) - vp_test;
			sout = sum(s1 .* s2, 2);
			%sout = max(0, sout);
			vars = sum((u' * y(ntrain,itest)).^2 + (v' * y(ntest,itest)).^2,2)/2;
			if k==length(ndims1)
				subplot(1,2,1),
				semilogx((results.cov_neur(1:npc)-sout)./vars)
				ypred = [u*vp_train; v*vp_test];
				vpred = vp_test;
				hold all;
				drawnow;
			end
			
			results.cov_res_beh_t(:,k,j,l) = gather_try(sout);
		end
	
	
	subplot(1,2,2),
	bx=cumsum(results.cov_neur(1:npc,td) - results.cov_res_beh_t(1:128,5,j,l));
	semilogx(ndims0,bx(npc,:)/sum(results.var_neur(1:npc,td)))
	disp([tlag(j) max(bx(npc,:)/sum(results.var_neur(1:npc,td)))])
	ylim([0 max(bx(npc,:)/sum(results.var_neur(1:npc,td)))]);
	drawnow;
	end
	disp([tlag(j) max(bx(npc,:)/sum(results.var_neur(1:npc,td))) nansum(results.cov_neur_t(:,j),1)/nansum(results.var_neur)]);
	
end


%%
cov_neur_t = results.cov_neur_t;
var_neur = results.var_neur;
cov_res_beh_t = results.cov_res_beh_t;

