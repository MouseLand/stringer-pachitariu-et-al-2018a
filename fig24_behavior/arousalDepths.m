function arousalDepths(dataroot, matroot)

dall=load(fullfile(dataroot, 'dbspont.mat'));

for d = [1:length(dall.db)]
    %%
    dat = load(fullfile(dataroot,sprintf('spont_%s_%s.mat',dall.db(d).mouse_name,dall.db(d).date)));
    if isfield(dat.stat, 'redcell')
        Ff = dat.Fsp(~logical([dat.stat(:).redcell]),:);
        cellpos{d} = dat.med(~logical([dat.stat(:).redcell]),:);
    else
        Ff = dat.Fsp;
        cellpos{d} = dat.med;
    end
    cellpos{d} = cellpos{d}(sum(Ff,2)>0,:);
    Ff = Ff(sum(Ff,2)>0,:);
    
    
    %%
    tdelay = 1;
    y    = Ff(:,(tdelay+1):end);
    [NN,NT] = size(y);
    disp('>>>>>>>>>>>>> ');
    disp([NN NT]);
    
    % 1.2 second bins
    if dat.db.nplanes == 10
        tbin = 4;
    else
        tbin = 3;
    end
   
    y = bin2d(y, tbin, 2);
    whisk = dat.beh.whisker.motionSVD(:,1);
    whisk = whisk * sign(corr(dat.beh.runSpeed(:), whisk(:)));
    
    x = [dat.beh.runSpeed(:) dat.beh.pupil.area(:) whisk];
    x    = x(1:end-tdelay,:);
    x = bin2d(x, tbin, 1);
    
    % correlate x and y
    cc{d} = corr(y',x);
    
    
end

save(fullfile(matroot,'arousalcorr.mat'),'cc','cellpos');















