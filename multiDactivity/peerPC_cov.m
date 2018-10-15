function peerPC_cov(dataroot,matroot,useGPU,dex)

nPC = 2.^[0:10];
ndim0 =  1:512;

dall=load(fullfile(dataroot, 'dbspont.mat'));

clf;
%hold all;

%%
cov_neur = [];
var_all_neur = [];
rng('default');
%%
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
    
    npc = 1024;%min(size(vtrain,2));
    
	cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
	[u,s,v] = svdecon(cov);
	s1 = u(:,1:npc)' * Ff(ntrain,itest);
	s2 = v(:,1:npc)' * Ff(ntest,itest);
	sneur = sum(s1 .* s2, 2);
	varneur = sum(s1.^2 + s2.^2,2)/2;
	semilogx(sneur./varneur)
		
	cov_neur(:,d) = gather_try(sneur);
	var_neur(:,d) = gather_try(varneur);
	drawnow;
	
	if d == dex
		exampleV1 = gather_try(s1);
		exampleV2 = gather_try(s2);
	end
end

%%

save(fullfile(matroot,'PCpredict.mat'),'cov_neur','var_neur','exampleV1','exampleV2');