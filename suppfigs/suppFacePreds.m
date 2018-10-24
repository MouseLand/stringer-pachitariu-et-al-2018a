close all;
default_figure([15 1 9 9]);

%%
tstart = [1000 190 1 1850 300 1 200 300 1];

clf;
dp = [1:9];
for d = 1:length(dp)
    hp=my_subplot(6,3,3*floor((d-1)/3) + d,[.9 .9]);
    hp(1).Position(2) = hp(1).Position(2) - .01; 
    trange = tstart(d)+[1:750];
    if dp(d)==2
        trange = [1:200 290+[1:550]];
    end
    spks = ytest{dp(d)}(:,trange);
    stdspks = max(1e-2,std(ytest{dp(d)}, 1, 2));
    spks = spks ./ stdspks;
    spks = spks(isortembed{dp(d)},:);
    spks = my_conv2(spks, 6, 1);
    imagesc(spks,[0 .3])
    hold all;
    plot([1 60/1.2],[0 0]+size(spks,1)+100,'k');
    plot(-10+[0 0],size(spks,1)+[-1000 0],'k');
    if d==1
        ht=text(-.05,0, '1000 neurons','fontsize',8,'fontangle','normal');
        ht.Rotation=90;
    end
    box off;
    axis tight;
    axis off;
    title(sprintf('recording %d',dp(d)),'fontweight','normal')
    
    hp=my_subplot(6,3,3*(floor((d-1)/3))+3 + d,[.9 .9]);
    hp.Position(2) = hp.Position(2) + 0.01;
    hp.Position(2) = hp.Position(2) - .01; 
    spks = ypred{dp(d)}(:,trange);
    spks = spks ./ stdspks;
    spks = spks(isortembed{dp(d)},:);
    spks = my_conv2(spks, 1, 1);
    imagesc(spks,[0 0.3]);
    hold all;
    plot([1 60/1.2],[0 0]+size(spks,1)+100,'k');
    plot(-10+[0 0],size(spks,1)+[-1000 0],'k');
    if d==1
        text(0,0,'1 min','fontsize',8,'fontangle','normal');
    end
    box off;
    axis tight;
    axis off;
    
end
cm=colormap('gray');
colormap(flipud(cm));

%%
print(fullfile(matroot,'suppFacePreds.pdf'),'-dpdf');