function fig_traces(sproj,tshared,tproj,tshared2,tproj2,cm,dobin,titles)

if dobin
    tshared = bin2d(tshared(:),4);
    sproj   = bin2d(sproj, 4, 2);
    tproj   = bin2d(tproj, 4, 2);
    
    tshared2 = bin2d(tshared2(:),4);
    tproj2   = bin2d(tproj2, 4, 2);
    
end

hold all;
for j = 1:size(sproj,1)
    T = sproj(j,:);
    %         T = my_conv2(T, 1, 2);
    plot(-j*4.5 + zscore(T)*1, 'linewidth', 1, 'color', cm(j,:),'linewidth',.5)
    
end



plot(zscore(tshared)*1.4-(j+2.5)*4.5, 'linewidth',1,'color',.5*[1 1 1],'linewidth',.5);
text(10,-(j+3)*4.5+10,{'stim-behav','shared dimension'},'fontsize',8,'color',.4*[1 1 1],...
    'units','data');

ny = size(sproj,1) + 2.5;
for j = 1:size(tproj,1)
    T = tproj(j,:);
    %         T = my_conv2(T, 1, 2);
    plot(-ny*4.5-(j)*4.5 + zscore(T)*.8,'color',[.4 .6 1], 'linewidth', .5)
    hold all
end
plot(zscore(tshared2)*1.4-(ny+j+2)*4.5, 'linewidth',1,'color',.5*[1 1 1],'linewidth',.5);
text(10,-(ny+j+2)*4.5+9,{'stim-spont','shared dimension'},'fontsize',8,'color',.4*[1 1 1],...
    'units','data');

ny = size(sproj,1) + size(tproj,1)+4.5;

for j = 1:size(tproj2,1)
    T = tproj2(j,:);
    %         T = my_conv2(T, 1, 2);
    plot(-ny*4.5-(j)*4.5 + zscore(T)*.8, 'k', 'linewidth', .5)
    hold all
end


ht=text(-.04,.8,titles{1},'horizontalalignment','center','fontsize',8,'color',[.8 .3 .8],...
    'fontangle','normal');
ht.Rotation = 90;

ht=text(-.04,.42,titles{2},'color',[.4 .6 1],'horizontalalignment','center','fontsize',8,'fontangle','normal');
ht.Rotation = 90;

ht=text(-.04,.1,titles{3},'horizontalalignment','center','fontsize',8,'fontangle','normal');
ht.Rotation = 90;


axis tight
box off
axis off;
