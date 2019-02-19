load(fullfile(matroot,'stimvar_with_decode.mat'));

%%
close all;
default_figure([1 1 3.25 3.25]);

%%
clf;


cs = [];
cs(1,:) = [.5 .2 .7];
cs(2,:) = [1 .3 1];
cs(3,:) = [.8 .3 .8];
cs(4,:) = [.5 0 .5];
cs(5,:) = [.7 .4 .7];

yh=.03;

ktype = 1; % behav prediction
tpts = .8e4 + .03e4 + [1:.15e4];%2.09e4;
i=1;
hs{i}=my_subplot(1,1,1,[.93 .8]);
hs{i}.Position(2) = hs{i}.Position(2)+.05;

sdiff = diff(stimtpt);
son = find(sdiff==1);
soff = find(sdiff==-1);

tshared1 = tshared{1} * sign(mean(Ushared{1,1}));
tshared2 = tshared{2} * sign(mean(Ushared{2,2}));

%hs{i}.Position(3) = .5;
dobin=1;
tproj1=tprojS{1}(1:3,:);
tproj1= sign(skewness(tproj1,1,2)) .* tproj1;
tproj2=tprojS{2}(1:3,:);
tproj2= sign(skewness(tproj2,1,2)) .* tproj2;
titles = {'stim-only','behav-only','spont-only'};
fig_traces(sprojS{ktype}(1:2:10,tpts),tshared1(tpts),tproj1(:,tpts),tshared2(tpts),tproj2(:,tpts),cs,dobin,titles);

axes('position',[hs{i}.Position(1) hs{i}.Position(2)-yh*1.7 hs{i}.Position(3) yh*1.5]);
fig_stimon(stimtpt(tpts),0,.8);

axes('position',[hs{i}.Position(1) hs{i}.Position(2)-5*yh hs{i}.Position(3) yh*2.8]);
plot(bin2d(runspeed(tpts),4),'color',[.2 .8 .2]);
text(.03,1.1,'running','color',[.2 .8 .2],'fontangle','normal','fontsize',6);
hold all;
plot( [0 5*5/1.3],-6 * [1 1],'k','linewidth',2);
text(0.02,-.03,'5 s','horizontalalignment','center','fontsize',6,'fontangle','normal');
axis off;
axis tight;
set(gca,'fontsize',6);

%%
print(fullfile(matroot, 'suppSSzoom.pdf'),'-dpdf');