clearvars -except dat

nPC = 2.^[0:11];
%%
% clear Npca Npairs
clf;
hold all;
lambda = 10; % regularizer

load('../dbspont.mat');
%%
rng('default');
for d = 1%[1:length(db)]
    %%
    dat = load(sprintf('../spont_%s_%s.mat',db(d).mouse_name,db(d).date));
    
    if isfield(dat.stat, 'redcell')
        Ff = dat.Fsp(~logical([dat.stat(:).redcell]),:);
        med = dat.med(~logical([dat.stat(:).redcell]),:);
    else
        Ff = dat.Fsp;
        med = dat.med;
    end
    
    % xyz distance
    pd = ((med(:,1)-med(:,1)').^2 + (med(:,2)-med(:,2)').^2 + (med(:,3)-med(:,3)')).^.5;
    
    reps = 100 * (~(pd<70));
    
    %% bin spikes in 1.2 second bins
    if dat.db.nplanes==10
        tbin = 4;
    else
        tbin=3;
    end
    [NN, NT] = size(Ff);
    disp(NN);
    Ff    = squeeze(mean(reshape(Ff(:,1:floor(NT/tbin)*tbin),...
        NN, tbin, []),2));
    Ff = (Ff - mean(Ff,2));
    NT = size(Ff,2);
    %%
    rng(11212);
    isp     = randperm(NN);
    isp_tr  = isp(1:ceil(NN/2));
    isp_tst = isp(1+ceil(NN/2):NN);
    
    Lblock = 60;
    fractrain = 0.75;
    [itrain, itest] = splitInterleaved(NT, Lblock, fractrain, 1);
    
    %if rotateFf==1
    %    [Urand, Sv] = svd(gpuArray.randn(numel(isp_tr), 'single'));
    %    Ff(:, isp_tr) = Ff(:, isp_tr) * gather(Urand);
    %end
    %       %%
    tic;
    Ff = gpuArray(single(Ff));
    nneur = 2.^[12];
    for nn = 1:length(nneur)
        for i = 1:1000%length(isp_tst)
            isp_valid   = isp_tr(reps(isp_tst(i), isp_tr)>0);
            isp_valid   = isp_valid(1:min(numel(isp_valid),nneur(nn)));
            X           = Ff(isp_valid, itrain) / sqrt(size(Ff,1)) * sqrt(1e4);
            [U Sv V]    = svdecon(X);
            V           = bsxfun(@times, V , diag(Sv)'); %Ff(itrain,isp_valid) * U;
            
            %%
            for k = 1:length(nPC(nPC<=size(V,2)))
                W = (Ff(isp_tst(i), itrain) * V(:,1:nPC(k)))/...
                    (V(:,1:nPC(k))'*V(:,1:nPC(k)) + ...
                    lambda * gpuArray.eye(nPC(k), 'single'));
                %
                Wf = U(:,1:nPC(k)) * W';
                Fpred = Wf' * Ff(isp_valid, itest);
                rez = mean((Ff(isp_tst(i), itest) - Fpred).^2, 2);
                tV  = mean(Ff(isp_tst(i), itest).^2, 2);
                
                results.resid{d}(i,k,nn) = gather(rez);
                results.totvar{d}(i,k,nn) =  gather(tV);
            end
            if rem(i,100)==0
                toc
            end
        end
        
        %%
        plot(1-nanmean(results.resid{d}(1:i,:,nn),1)./nanmean(results.totvar{d}(1:i,:,nn)),'k')
        title(max(1-nanmean(results.resid{d}(1:i,:,nn),1)./nanmean(results.totvar{d}(1:i,:,nn))));
        hold all;
        drawnow;
    end
    
end

%%
save('peerPCA/PCApred2.mat','results','expv_neurons');
%
% save('results/PCApred.mat', 'eVarTst', 'eVarTra', 'TotVarF', 'MaxVarPerc')
%%
expv_neurons=[];

for d = [1:numel(results.resid)]
    expv_neurons = cat(2, expv_neurons, (1-nanmean(results.resid{d}(:,1:11))./nanmean(results.totvar{d}(:,1:11)))');
    %expv
    
end


