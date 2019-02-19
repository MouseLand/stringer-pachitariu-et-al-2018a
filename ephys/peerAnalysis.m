function [expv, resid, totvar, nPC] = peerAnalysis(stall, iprobe, Wh, dmin)

lambda = 10; % regularizer

% exclude neurons on the same probe within x um
pd = max(1,abs(Wh - Wh'));

pd = pd .* (max(1,abs(iprobe - iprobe')*(dmin+1)));

reps = ~(pd<dmin);
    
Ff = stall - mean(stall,2);
[NN,NT] = size(Ff);
    
%%
rng(11212);
isp     = randperm(NN);
isp_tr  = isp(1:ceil(NN/2));
isp_tst = isp(1+ceil(NN/2):NN);

Lblock = 60;
fractrain = 0.75;
[itrain, itest] = splitInterleaved(NT, Lblock, fractrain,1);

nPC = 2.^[0:9];
%if rotateFf==1
%    [Urand, Sv] = svd(gpuArray.randn(numel(isp_tr), 'single'));
%    Ff(:, isp_tr) = Ff(:, isp_tr) * gather(Urand);
%end
%       %%
tic;
Ff = gpuArray(single(Ff));
resid = NaN*zeros(NN, length(nPC),'single');
totvar = NaN*zeros(NN, length(nPC),'single');
for i = 1:length(isp_tst)
    isp_valid   = isp_tr(reps(isp_tst(i), isp_tr)>0);
    isp_valid   = isp_valid(1:min(numel(isp_valid),2000));
    X           = Ff(isp_valid, itrain);
    [U Sv V]    = svdecon(X);
    V           = bsxfun(@times, V , diag(Sv)');
    
    %%
    for k = 1:length(nPC)
        W = (Ff(isp_tst(i), itrain) * V(:,1:nPC(k)))/...
            (V(:,1:nPC(k))'*V(:,1:nPC(k)) + ...
            lambda * gpuArray.eye(nPC(k), 'single'));
        %
        Wf = U(:,1:nPC(k)) * W';
        Fpred = Wf' * Ff(isp_valid, itest);
        rez = mean((Ff(isp_tst(i), itest) - Fpred).^2, 2);
        tV  = mean(Ff(isp_tst(i), itest).^2, 2);
        
        resid(isp_tst(i),k) = gather(rez);
        totvar(isp_tst(i),k) =  gather(tV);
    end
    if rem(i,100)==0
        toc
    end
end
    
expv = 1 - nanmean(resid,1) ./ nanmean(totvar,1);

%%
clf;
plot(expv,'r')
title(max(expv));
drawnow;



