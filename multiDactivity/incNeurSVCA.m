function incNeurSVCA(dataroot,matroot,useGPU,dex)

dall=load(fullfile(dataroot, 'dbspont.mat'));

nneur0 = 2.^[8:12];
npc = 1024;

clf;
%%

cov_neur = NaN*zeros(npc, length(nneur0), 10, length(dall.db));
var_neur = NaN*zeros(npc, length(nneur0), 10, length(dall.db));
nneurs = [];

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
    
	nneur1 = nneur0;
	nmax = min(length(ntest),length(ntrain));
	if nmax < nneur1(end)
		nneur1(end) = nmax;
	end
	
	for n = 1:length(nneur1)
		% do multiple subsets and avg
		for it = 1:10
			ntrain1 = ntrain(randperm(length(ntrain), nneur1(n)));
			ntest1 = ntest(randperm(length(ntest), nneur1(n)));
			[sneur, varneur, u, v] = SVCA(Ff, min(npc, nneur1(n)), ntrain1, ntest1,	itrain, itest);
			cov_neur(1:length(sneur),n,it,d) = gather_try(sneur);
			var_neur(1:length(sneur),n,it,d) = gather_try(varneur);
		
		end
		semilogx(cov_neur(:,n,it,d)./var_neur(:,n,it,d))
		
		ylim([0 1]);
		hold all;
		drawnow;
		
	end
	nneurs(:,d) = nneur1;
end

%%

save(fullfile(matroot,'incNeurSVCA.mat'),'cov_neur','var_neur','nneurs');