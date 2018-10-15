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
    
    [utrain,~,vtrain1] = svdecon(Ff(ntrain,itrain));
    [utest,~,vtest1] = svdecon(Ff(ntest,itrain));
    npc = 1024;%min(size(vtrain,2));
    vtrain1 = vtrain1(:,1:npc);
    vtest1 = vtest1(:,1:npc);
    utrain = utrain(:,1:npc);
    utest = utest(:,1:npc);
	vtrain2 = Ff(ntrain,itest)' * utrain;
    vtest2 = Ff(ntest,itest)' * utest;
    %%
	%cov = vtrain1' * vtest1;
	%[u,s,v] = svdecon(cov);
	%sout = u' * vtrain2' * vtest2 * v;
	cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
	[u,s,v] = svdecon(cov);
	s1 = u' * Ff(ntrain,itest);
	s2 = v' * Ff(ntest,itest);
	sneur = sum(s1 .* s2, 2);
	varneur = sum(s1.^2 + s2.^2,2)/2;
	semilogx(sneur./varneur)
		
	cov_neur(:,d) = sneur(1:1024);
	var_all_neur(:,d) = varneur(1:1024);
	drawnow;
end

%%

save(fullfile(matroot,'PCpredict.mat'),'expvPC','expvPCCov','expvPCVar','expv128','ndims1','exampleV1','exampleV2');