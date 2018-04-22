function peerExcludeNeighbors(dataroot,matroot,useGPU)

nPC = 2.^[0:11];

dall=load(fullfile(dataroot, 'dbspont.mat'));

clf;
hold all;
lambda = 10; % regularizer

%%
rng('default');
for d = [1:length(db)]
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
    
    tic;
    if useGPU
        Ff = gpuArray(single(Ff));
    end
    nneur = 2.^[12];
    for nn = 1:length(nneur)
        for i = 1:min(numel(isp_tst), 1000)
            isp_valid   = isp_tr(reps(isp_tst(i), isp_tr)>0);
            isp_valid   = isp_valid(1:min(numel(isp_valid),nneur(nn)));
            X           = Ff(isp_valid, itrain) / sqrt(size(Ff,1)) * sqrt(1e4);
            [U Sv V]    = svdecon(X);
            V           = bsxfun(@times, V , diag(Sv)');
            
            %%
            for k = 1:length(nPC(nPC<=size(V,2)))
                if useGPU
                    W = (Ff(isp_tst(i), itrain) * V(:,1:nPC(k)))/...
                        (V(:,1:nPC(k))'*V(:,1:nPC(k)) + ...
                        lambda * gpuArray.eye(nPC(k), 'single'));
                else
                    W = (Ff(isp_tst(i), itrain) * V(:,1:nPC(k)))/...
                        (V(:,1:nPC(k))'*V(:,1:nPC(k)) + ...
                        lambda * eye(nPC(k), 'single'));
                end
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

expv_neurons=[];

for d = [1:numel(results.resid)]
    expv_neurons = cat(2, expv_neurons, (1-nanmean(results.resid{d}(:,1:11))./nanmean(results.totvar{d}(:,1:11)))');
end

save(fullfile(matroot,'PCApred.mat'),'results','expv_neurons');



