function suppDists(matroot)

pc = load(fullfile(matroot,'corr1stpc.mat'));
clust = load(fullfile(matroot,'clust1D.mat'));

load(fullfile(matroot,'expv_behavior_neurons.mat'));

%%

figX = 4;
figY = 7;
close all;
default_figure([1 1 figX figY]);


%%
clf;
nC = size(clust.results.depth,1);
cg = colormap('jet');
rng(6);
cg = cg(randperm(64),:);
cg = cg(round(linspace(1,64,nC)),:);

ndat = numel(pc.results.runcorr);
load('cdat.mat');

xh = 0.6;
yh = 0.6;

dex = 2;

i=0;
% --------- correlation vs distance
i = i+1;
hs{i} = my_subplot(4,2,1,[xh yh]);
dx = pc.results.dbins(1:end-1) + diff(pc.results.dbins(1:2));
hold all;
for j = [1:size(pc.results.dmean,2)]
    plot(dx, pc.results.dmean(:,j),'color',cdat(j,:),'linewidth',.25)
end
j=dex;
plot(dx, pc.results.dmean(:,j),'color',cdat(j,:),'linewidth',1)
xlabel('distance (μm)','fontname','Arial','fontsize',7)
ylabel('mean correlation');
ylim([-.1 .1]);
axis square;

i=i+1;
hs{i}=my_subplot(4,2,2,[xh yh]);
axis off;
hp=hs{i}.Position;
axes('Position',[hp(1)-.01 hp(2) hp(3:4)*1.12]);
hold all;
iplane = 5;
for j = 1:nC
    plot(clust.results.pos(clust.results.iclust==j & clust.results.pos(:,3)==iplane*35,1),...
        clust.results.pos(clust.results.iclust==j & clust.results.pos(:,3)==iplane*35,2),'.',...
        'color',cg(j,:),'markersize',4);
end
plot([0 500],[0 0]-30,'k','linewidth',2)
text(0,0,'500 \mum','fontsize',8,'fontangle','normal');
axis tight;
box off;
axis square;
axis off;
%text(-.24,1.35, 'example plane','fontangle','normal');

i=i+1;
hs{i}=my_subplot(4,2,3,[xh yh]);
hold all;
for j = 1:nC
    plot([clust.results.indist(j,dex) clust.results.outdist(j,dex)],'color',cg(j,:));
end
axis([.75 2.25 0 550])
box off;
axis square;
ylabel('pairwise distance (μm)','fontname','Arial','fontsize',7)
set(gca,'xtick',[1 2],'xticklabel',{'in-group','out-of-group'});

i=i+1;
hs{i}=my_subplot(4,2,4,[xh yh]);
hold all;
for d = 1:ndat
    plot(mean([clust.results.indist(:,d) clust.results.outdist(:,d)],1),'color',cdat(d,:),'linewidth',.5);
end
axis([.75 2.25 0 550])
box off;
axis square;
ylabel('pairwise distance (μm)','fontname','Arial','fontsize',7)
set(gca,'xtick',[1 2],'xticklabel',{'in-group','out-of-group'});


% cm is running pupil whisker
cm(1,:) = [.2 .8 .2];
cm(2,:) = [.5 .6 .5];
cm(3,:) = [0 .2 0];

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
ca = repmat(cm,2,1);

i=i+1;
hs{i}=my_subplot(4,2,5,[xh yh]);
hold all;
for j = 1:nC
    plot([clust.results.indepth(j,dex) clust.results.outdepth(j,dex)],'color',cg(j,:));
end
axis([.75 2.25 0 150])
box off;
axis square;
ylabel('distance in depth (μm)','fontname','Arial','fontsize',7)
set(gca,'xtick',[1 2],'xticklabel',{'in-group','out-of-group'});

i=i+1;
hs{i}=my_subplot(4,2,6,[xh yh]);
hold all;
for d = 1:ndat
    plot(mean([clust.results.indepth(:,d) clust.results.outdepth(:,d)],1),'color',cdat(d,:),'linewidth',.5);
end
axis([.75 2.25 0 150])
box off;
axis square;
ylabel('distance in depth (μm)','fontname','Arial','fontsize',7)
set(gca,'xtick',[1 2],'xticklabel',{'in-group','out-of-group'});

i=i+1;
hs{i}=my_subplot(4,2,7,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)+0.015;
axis([0.5 7 deps(1)+150 deps(end)+150]);
ylabel('depth (μm)','fontname','Arial','fontsize',7)
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
        ha=area(bd,nb,'edgecolor','none','facecolor',ca(j,:));
    else
        ha=area(bd,nb,'edgecolor','none','facecolor',min(1,ca(j,:)+.2));
    end
    ha.ShowBaseLine = 'off';
    view([90 90]);
    axis off;
    axis tight;
    xlim([deps(1) deps(end)]);
   
end


tstr{1} = 'correlation half 1';
tstr{2} = 'correlation half 2';
tstr{3} = '';
tstr{4} = '';

% -------------- LETTERS
hp=.12;
hy=1.25;
deffont=8;
for j = [1:length(hs)]
    hp0=hp;
    hy0 = hy;
    
    hpos = hs{j}.Position;
    lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
    axes('position', lpos);
    text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
        
    if j==10
        tsh=.02;
    else
        tsh=.03;
    end
    axes('position', [lpos(1)+tsh lpos(2)-.006 lpos(3) lpos(4)]);
    %text(0,0, tstr{j},'fontangle','normal','fontsize',8);
    axis([0 1 0 1]);
    axis off;
    
end


%%
print(fullfile(matroot,'suppDists.pdf'),'-dpdf')

% 
% % PC DISTANCE
% % ----------- example plane -------%
% i = i+1;
% hs{i} = my_subplot(3,3,2,[xh yh]);
% axis off;
% hp=hs{i}.Position;
% %hp(1)=hp(1)-.045;
% %hp(2)=hp(2)-.015;
% %hp(3)=hp(3)+.035;
% %hp(4)=hp(4)+.035;
% axes('position',hp);
% iplane = 35*2;
% icell = find(pc.results.cellmed(:,3)==iplane);
% med = pc.results.cellmed(icell,1:2);
% ctype= pc.results.cellpc(icell);
% hold all;
% cm=[1 0 0; 0 0 1];
% NN = size(med,1);
% for j = randperm(NN)
%     plot(med(j,1),med(j,2),'.','markersize',3,'color',cm(ctype(j),:));
% end
% axis off;
% axis square;
% plot([0 500],[0 0]-30,'k','linewidth',2)
% text(0,0,'500 \mum','fontsize',8,'fontangle','normal');
% text(.45, 0, {'positive weight'},'color',cm(1,:),'fontsize',8,'fontangle','normal');
% text(.45, -.1, {'negative weight'},'color',cm(2,:),'fontsize',8,'fontangle','normal');
% 
% i = i+1;
% hs{i} = my_subplot(3,3,3,[xh yh]);
% hold all;
% for d = 1:ndat
%     lw = .5;
%     if d==pc.dex
%         lw = 1;
%     end
%     plot([pc.results.indist(:,d) pc.results.outdist(:,d)]','color',cdat(d,:),'linewidth',lw)
% end
% p = signrank(pc.results.indist(:),pc.results.outdist(:));
% ylim([0 580]);
% xlim([.75 2.25]);
% box off;
% ylabel('pairwise distance (um)');
% set(gca,'xtick',[1 2],'xticklabel',{'same sign',['opp sign']});
% axis square;
