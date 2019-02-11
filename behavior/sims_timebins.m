
x = randn(2e4,500);
x = x';
x = my_conv2(x, 60, 2);
nPC = size(x,1);

NN = 8e3;
[a,~,b] = svdecon(exprnd(1,NN,nPC));
ncomps=32;
a = a(:,1:ncomps);
s = diag(1./[1:32].^.5);
b = b(:,1:ncomps);

z = b' * x;
z = zscore(z, 1, 2);
y = a * s * z;

%% add noise
%y = (y + .02 * randn(size(y))) .* max(0, 0.1 + 5*exprnd(1,size(y)));
%y = max(0, y / std(y(:)) + 10 * randn(size(y)));
%y = y .* max(0, 0.1 + 10*exprnd(1,size(y)));
%y = y + exprnd(3,size(y));

%
%yp = y;

yp = poissrnd(max(0, y + 0.02*randn(size(y))) * 2);
% firing rate
mean(yp(:))/0.03

clf;
hold all;
plot(x(1,:));
plot(my_conv2(yp(1,:),5,2))

%%

%
tbin = 10;
Ff  = double(bin2d(yp, tbin, 2));
xb  = double(bin2d(x, tbin, 2));

Ff = Ff - mean(Ff,2);
NT = size(Ff,2);
ntrain = 1:NN/2;
ntest = NN/2+1:NN;

% divide time in half
%Ff = randn(size(Ff));
Lblock = 60;
fractrain = 0.5;
[itrain, itest] = splitInterleaved(NT, Lblock, fractrain, 1);
tic;

[sneur, varneur, u, v] = SVCA(Ff, 1024, ntrain, ntest, itrain, itest);
npc = 128;
u = u(:,1:npc);
v = v(:,1:npc);

clf;
%subplot(1,2,1);
semilogx(sneur./varneur)
td = 1;
d = 1;
results.cov_neur(1:length(sneur),td,d) = gather_try(sneur);
results.var_neur(1:length(sneur),td,d) = gather_try(varneur);
%
%%
lam = .01;

clf;
ndims1 = 2.^[0:6];
[atrain, btrain] = CanonCor2(Ff(ntrain,itrain)'*u, xb(:,itrain)', lam);
[atest, btest]   = CanonCor2(Ff(ntest,itrain)'*v, xb(:,itrain)', lam);
k=0;
for n = ndims1
	k=k+1;
	vp_train     = atrain(:,1:n) * btrain(:,1:n)' * xb(:,itest);
	vp_test      = atest(:,1:n) * btest(:,1:n)' * xb(:,itest);
	
	s1 = u' * Ff(ntrain,itest) - vp_train;
	s2 = v' * Ff(ntest,itest) - vp_test;
	sout = sum(s1 .* s2, 2);
	vars = sum((u' * Ff(ntrain,itest)).^2 + (v' * Ff(ntest,itest)).^2,2)/2;

	semilogx((results.cov_neur(1:npc,td,d)-sout)./vars)
	hold all;
	
	ypred = [u*vp_train; v*vp_test];
	pred = 1 - nanmean(nanmean((ypred - Ff(:,itest)).^2,2))/nanmean(nanmean(Ff(:,itest).^2));
	neurpred = 1 - nanmean((ypred - Ff(:,itest)).^2,2)./nanmean(Ff(:,itest).^2,2);
	fprintf('face pred: %2.2f\n',pred);
	vpred = vp_test;
	
	results.cov_res_beh(:,k,td,d) = max(0,gather_try(sout));
	results.var_beh(:,k,td,d) = gather_try(vars);
	ylim([0 1]);
	xlim([0 128]);
	drawnow;
	
end

expcbeh = cumsum(results.cov_neur(1:npc,td,d) - results.cov_res_beh(:,:,td,d))./cumsum(results.var_neur(1:npc,td,d));
expcbeh = squeeze(expcbeh(128,:));
expcbeh


%%
ops.nCall = [30 100]; % number of clusters
ops.iPC = 1:200; % PCs to use
ops.useGPU = useGPU; % whether to use GPU
ops.upsamp = 100; % upsampling factor for the embedding position
ops.sigUp = 1; % stddev for upsampling

[isort, ~, Sm] = mapTmap(Ff, ops);
			
%%
clf;
subplot(2,1,1),
it = find(itest)
imagesc(zscore(my_conv2(Ff(isort,it(1:1000)),[10 2],[1 2]),1,2),[-1 3])

subplot(2,1,2),
imagesc(zscore(my_conv2(ypred(isort,1:1000),[10 2],[1 2]),1,2),[-1 3])
