function fig2(matroot)


load(fullfile(matroot,'example_behavior.mat'));
load(fullfile(matroot,'expv_behavior_neurons.mat'));

%%
default_figure([1 1 6 8]);
%%
tstr   = {'running','pupil area', 'whisking', 'running+pupil','running+whisker',...
    'pupil+whisker','running+pupil+whisker'};


% cm is running pupil whisker
cm(1,:) = [.2 .8 .2];
cm(2,:) = [.5 .6 .5];
cm(3,:) = [0 .2 0];

cpc = max(0,colormap('spring'));
cpc = cpc(1:32,:);
nPCs = 4;
cpc = cpc(round(linspace(1,size(cpc,1),nPCs)),:);
ndat = size(expv_behavior,2);
load('cdat.mat');

% independent contributions
indvar = squeeze(expv_behavior(3,:,7) - expv_behavior(2,:,4:6));
indvar = indvar(:, [3 2 1]);

i = 0;

clf;
trange = 7300+[1:500];
i=i+1;
hp=my_subplot(4,3,1,[.6 .6]);
axis off;
hs{i} = my_subplot(4,1,1,[.87 .55]);
hs{i}.Position(1)=hp.Position(1);
%hs{i}.Position(3)=hp.Position(3);
xp = zscore(xex(trange,:),1,1);
xp = xp - min(xp,[],1);
xp(:,1) = xp(:,1)/1.3;
xp = min(5, xp);
%xp = xp ./ max(xp,[],1);
ht=text(-.07,0.07,'normalized units','fontsize',8,'fontweight','bold','fontangle','normal');
ht.Rotation=90;
ht=text(-.07,-1.55,'correlation','fontsize',8,'fontweight','bold','fontangle','normal');
ht.Rotation=90;
%ht=text(-.06,-5.1,'weights','fontsize',8,'fontweight','bold','fontangle','normal');
%ht.Rotation=90;
for j = [3 1 2]
    plot([0:1/3:1/3*(numel(trange)-1)],xp(:,j),'color',cm(j,:),'linewidth',1);
    hold all;
    text(.02,1.-(j-1)*.12,tstr{j},'color',cm(j,:),'fontweight','bold',...
        'fontsize',8,'fontangle','normal');
end
xlabel('time (s)');
axis tight;
box off;
set(gca,'color','none','YColor','none');
%ylabel('normalized units');

i=i+1;
hs{i} = my_subplot(4,3,4,[.6 .6]);
hs{i}.Position(2) = hs{i}.Position(2) + .02;
cmat = mean(cc,3);
cmat = cmat - diag(NaN*diag(cmat));
cmat = [squeeze(cc(1,2,:)) squeeze(cc(1,3,:)) squeeze(cc(2,3,:))];
for j = 1:ndat-1
    plot(cmat(j,:),'color',cdat(j,:))
    hold all;
end
plot(mean(cmat(1:end-1,:),1),'ko-','linewidth',2,'markerfacecolor','k','markersize',4)
axis tight;
box off;
ylim([0 .85]);
set(gca,'xtick',[1 2 3],'xticklabel',...
    {'run-pupil','run-whisk','pupil-whisk'});
ylabel('correlation');
axis square;
xtickangle(45);
xlim([0.75 3.25]);

% pcolor([cmat nan(3,1);nan(1,4)])
% shading flat;
% set(gca,'ydir','reverse');
% colormap(hs{i},'jet');
% set(gca,'xtick',[1 2 3]+.5,'xticklabel',{'running','pupil','whisking'});
% set(gca,'ytick',[1 2 3]+.5,'yticklabel',{'running','pupil','whisking'});
% colorbar;
% title('correlation');

i=i+1;
hs{i} = my_subplot(4,3,5,[.6 .6]);
hs{i}.Position(2) = hs{i}.Position(2) + .02;
%hs{i}.Position(3)=hs{i}.Position(3)+.03;
ev=[squeeze(expv_behavior(1,:,1:3)) squeeze(expv_behavior(3,:,7))'];
for j = 1:ndat
    plot([1:4],ev(j,:),'color',cdat(j,:))
    hold all;
end
plot([1:4],mean(ev),'ko-','linewidth',2,'markerfacecolor','k','markersize',4)
box off;
set(gca,'xtick',[1:8],'xticklabel',{'running','pupil area','whisking','all 3'});
xtickangle(45);
axis square;
axis tight;
ylabel({'variance explained','(test set)'});
xlim([0.5 4.5]);
ylim([0 0.065]);

i=i+1;
hs{i} = my_subplot(4,3,6,[.6 .6]);
hs{i}.Position(2) = hs{i}.Position(2) + .02;
for j = 1:ndat
    plot([1:3],indvar(j,:)','color',cdat(j,:))
    hold all;
end
plot([1:3],mean(indvar,1)','ko-','linewidth',2,'markerfacecolor','k','markersize',4)
set(gca,'xtick',[1:8],'xticklabel',{'running','pupil area','whisking','all 3'});
xtickangle(45);
axis square;
ylabel({'variance explained','(unique, test set)'});
box off;
ylim([0 .065]);
xlim([0.75 3.25]);
%

% --------------- SCHEMATIC
i=i+1;
hs{i} = my_subplot(4,1,3);
hs{i}.Position(1)=hs{1}.Position(1)-.02;
hs{i}.Position(3)=hs{1}.Position(3)*.25;
hs{i}.Position(2)=hs{i}.Position(2)+.02;
hs{i}.Position(4)=hs{i}.Position(4)*1.;

dex=2;
trange = 200+[1:100];
x = (xtest1d{dex}(:,trange) - nanmean(xtest1d{dex}(:,trange),2))...
    ./ nanstd(xtest1d{dex}(:,trange),1,2);
x = real(x);
%x(3,:) = x(3,:) * -1;
NT = numel(trange);
for j = 1:3
    plot(x(j,:)*.55-j*4.2-.5,'color',cm(j,:));
    hold all;
    slope = (-8 +j*4)/NT;
    text(0,.88-(j-1)/4,tstr{j},'color',cm(j,:),'fontsize',8);%,'fontweight','bold');
    %ht=text(2*NT, (-j*4 - 8)/2+1,'>','units','data',...
    %    'fontsize',12,'fontweight','bold');
    %ht.Rotation = slope*1.5*360;
    
end
plot([0 30/1.2],[1 1]*(-j*4.2-2),'k');
text(0,0.08,' 30 s','fontsize',8,'fontangle','normal');
axis off;
cp = 0*[1 1 1];
axis off;
axis tight;
ylim([-j*5-1 1]);
text(0,1.07,'behaviors',...
    'fontsize',8);


xpred0 = real((c1d(:,dex)' * real(xtest1d{dex}(:,:))));

hp=hs{i}.Position;
hp(1)=hp(1)+.28;
%hp(3)=.12;
hp(4)=.03;
axes('position',hp);
plot(xpred0(trange),'color',cp);
text(.5, 5.3,{'best 1D','weighted','average'},...
    'fontsize',8,'horizontalalignment','center');
axis tight;
axis off;

hp=hs{i}.Position;
hp(1)=hp(1)+.55;
hp(2) = hp(2)-.0;
hp(4) = hp(4)*.9;
axes('position',hp);
n1d = real(n1d);
nsc = max(n1d(:,dex))*2;
nd = n1d(:,dex)/nsc;
xsc =  nanstd(xpred0(trange));
xpred = xpred0(trange) / xsc;
hold all;
for j =1:nPCs
    vplot = double(vtest{dex}(j,trange)/nsc/xsc);
    plot(vplot - j*3.7+.5,'color',cpc(j,:))
    plot(nd(j) * xpred - j*3.7+.5,'color',cp);
    r = 1 - (var(vtest{dex}(j,:)' - n1d(j,dex)*xpred0') ...
        / var(vtest{dex}(j,:)));
    text(70,-j*3.7+3.6, sprintf('r^2=%1.2f',r),'color',cpc(j,:),...
        'units','data','fontweight','bold','fontsize',8);
end
text(0,1.15,'neural PCs',...
    'fontsize',8);

ht=text(1.4,1.15,'neurons',...
    'fontsize',8,'HorizontalAlignment','center');


axis tight;
axis off;

% --- BALLS
nHidden = 1;
call{1} = cm;
call{2} = cp;
call{3} = cpc;
pos(1) = hs{i}.Position(1)+hs{i}.Position(3)*1.2;
pos(3) = hs{i}.Position(3)*1.2;
pos(2) = hs{i}.Position(2)+.0;
pos(4) = hs{i}.Position(4)*.8;
nls = [3 1 nPCs];
msize = 14;
layeredNet(pos, nls, call, msize, [0 0 0]);
% little ellipses
axes('position',[pos(1) pos(2)-.044 pos(3) .02]);
ipos = [1 2 3];
hold all;
for j = 3
    plot(ones(1,3)*ipos(j),[1:3],'k.');
end
axis tight;
xlim([1 3]);
axis off;

%pos = hs{i}.Position;
xb0 = pos(1)+pos(3)*2;
yb0 = pos(2)+.002;
axes('position',[xb0 yb0 .09 pos(4)]);
ny = 10;
hold all;
for k = 1:ny
    xb=(j-1);
    yb=(k-1);
    %    if j==1
    for m = 1:nPCs
        yp = (ny-1)/3 * (m-1);
        plot([-4 xb],[yp yb],'k','linewidth',1);
        hold all;
    end
    %text(0.,3,'neurons','fontsize',12);
    %   end
    plot(xb,yb,'ko','markerfacecolor','w','linewidth',2,'color','k',...
        'markersize',7);
    
end
axis tight;
axis off;


% weights of best 1D
i=i+1;
hs{i}=my_subplot(4,3,10,[.6 .6]);
cnorm = normc(c1d);
for j = 1:ndat
    plot(cnorm(:,j),'color',cdat(j,:));
    hold all;
end
plot(mean(cnorm,2),'ko-','linewidth',2,'markerfacecolor','k','markersize',4)
plot([1 3],[0 0],'k','linewidth',0.5)
set(gca,'xtick',[1:8],'xticklabel',{'running','pupil area','whisking'});
axis square;
box off;
ylabel('weights');
xtickangle(45);
xlim([0.75 3.25]);

% 1D VS 3D
i=i+1;
hs{i}=my_subplot(4,3,11,[.6 .6]);
ev=[squeeze(expv_behavior(1,:,7))' squeeze(expv_behavior(3,:,7))'];
plot([0 .07],[0 .07],'k','linewidth',1);
hold all;
for j = 1:ndat
    plot(ev(j,1),ev(j,2),'wo','markerfacecolor',cdat(j,:))
    hold all;
end
box off;
title({'variance explained','(test set)'},'fontsize',8,'fontweight','normal');
xlabel('best 1D model');
ylabel('full 3D model');
axis([0 .07 0 .07]);
set(gca,'xtick',[0:.02:.06],'ytick',[0:.02:.06]);
grid on;
grid minor;
grid minor;


% plot of arousal with depth
ccall=[];
posall=[];
for d = 1:9
    ccall = [ccall; ccarousal{d}];
    posall = [posall; cellpos{d}];
end
%ic=max(abs(ccall),[],2)>.2;
%ccall = ccall(ic,:);
%posall = posall(ic,:);
cc2=[ccall ccall*-1];
[~,ix] = max(cc2,[],2);
ng = size(cc2,2);
hold all;
deps = sort(unique(posall(:,3)));
deps = [deps; deps(end)+35] - 17;
cg = repmat(cm,2,1);

i=i+1;
hs{i}=my_subplot(4,3,12,[.6 .6]);
axis([0.5 7 deps(1)+150 deps(end)+150]);
ylabel('depth (\mum)');
set(gca,'Ydir','reverse')
set(gca,'xtick',[1:6],'xticklabel',{'running +','pupil area +','whisking +',...
    'running -','pupil area -','whisking -'});
xtickangle(45);
grid on;
hp=hs{i}.Position;
hp0 = hp(3);
hp(3) = hp0/(ng+3);
hp(1) = hp(1) - hp0/(ng+.5)/2;% - hp0/ng + .005;
nball = histcounts(posall(:,3), deps);
for j = 1:ng
    hp(1) = hp(1) + hp0/(ng+.5);
    axes('position',hp);
    nb = histcounts(posall(ix==j,3), deps);
    nb = nb/sum(nb)./nball;
    bd = deps(1:end-1)+17;
    if j<=ng/2
        ha=area(bd,nb,'edgecolor','none','facecolor',cg(j,:));
    else
        ha=area(bd,nb,'edgecolor','none','facecolor',min(1,cg(j,:)+.2));
    end
    ha.ShowBaseLine = 'off';
    view([90 90]);
    axis off;
    axis tight;
    xlim([deps(1) deps(end)]);
    
end


%
% -------------- LETTERS
hp=.08;
hy=1.23;
deffont=8;
for j = [1:length(hs)]
    if j==5
        hp0 = hp-.02;
    else
        hp0=hp;
    end
    hpos = hs{j}.Position;
    axes('position', [hpos(1)-hp0 hpos(2)+hpos(4)*hy(1) .01 .01]);
    text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
end