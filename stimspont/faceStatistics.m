function faceStatistics(dataroot, matroot)

load(fullfile(dataroot,'dbstimspont.mat'));

vface=[];
for d = 1:length(dbs)
    db = dbs(d);
    dat = load(fullfile(dataroot,...
        sprintf('stimspont_%s_%s.mat',db.mouse_name,db.date)));
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
        vface(:,j,d) = (var(x{j},1,1))';
    end
    vface(:,3,d) = vtot(:);
    
    if d == 1
        smot = squeeze(sign(mean(mean(dat.beh.face.motionMask,1),2)));
        faceEx = dat.beh.face.motionSVD .* smot';
        stimtpt= logical(dat.stimtpt);
    end
    
end
    
%%
save(fullfile(matroot,'faceSpectrum.mat'),'vface','faceEx','stimtpt');
