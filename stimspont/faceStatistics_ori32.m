function faceStatistics_ori32(dataroot, matroot)

%%
load(fullfile(dataroot,'dbori32.mat'));
% first orientation dataset doesn't have spont periods
dbs = dbs(2:end);

%%
vface=NaN*ones(500,3,length(dbs));
for d = 1:length(dbs)
    db = dbs(d);
    dat = load(fullfile(dataroot,...
        sprintf('orispont_%s_%s.mat',db.mouse_name,db.date)));
    tnoface = isnan(dat.beh.face.motionSVD(:,1));
   
    fprintf('recording %d\n',d);
    
    % fit to periods with or without stim
    tstim = dat.stimtpt & ~tnoface;
    tspont  = ~dat.stimtpt & ~tnoface;

    clear x;
    x{1} = dat.beh.face.motionSVD(tspont,:);
    x{2} = dat.beh.face.motionSVD(tstim,:);
    
    vtot = var(dat.beh.face.motionSVD(~tnoface,:),1,1);
    
    for j = 1:2
        vface(1:size(x{1},2),j,d) = (var(x{j},1,1))';
    end
    vface(1:size(x{1},2),3,d) = vtot(:);    
end
    
%%
save(fullfile(matroot,'faceSpectrum_ori32.mat'),'vface');
