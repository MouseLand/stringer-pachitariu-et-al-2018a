function quantifyBehavior(dataroot, matroot, dex)

dall=load(fullfile(dataroot, 'dbspont.mat'));

cc = [];

%%
for d = [1:length(dall.db)]

    dat = load(fullfile(dataroot,sprintf('spont_%s_%s.mat',dall.db(d).mouse_name,dall.db(d).date)));
   
    whisk = dat.beh.whisker.motionSVD(:,1);
    whisk = whisk * sign(corr(dat.beh.runSpeed(:), whisk(:)));
    
    x = [dat.beh.runSpeed dat.beh.pupil.area(:) whisk];
    
    cc(:,:,d) = corr(x);
    
    if d == dex
        xex = x;
    end
    
end
cc = real(cc);

%%

save(fullfile(matroot,'example_behavior.mat'),'cc','xex');