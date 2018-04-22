load('../dbspont.mat');

dex = 2;
cc = [];

%%
for d = [1:length(db)]

    dat = load(sprintf('../spont_%s_%s.mat',db(d).mouse_name,db(d).date));
   
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

save('example_behavior.mat','cc','xex');