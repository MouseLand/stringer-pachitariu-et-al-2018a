function fig6(matroot)

load(fullfile(matroot,'faceSpectrum.mat'));
load(fullfile(matroot,'stimvar.mat'));
load(fullfile(matroot,'stimfaceRepresentation.mat'));
load(fullfile(matroot,'exResponses.mat'));

%%
close all;
default_figure([15 1 9 8.5]);

%%
rng('default');
kt = 1; % face prediction
tpts = .645e4 + [1:.92e4];%2.09e4;

dx = .25;
dy = .5;

sdiff = diff(stimtpt);
son = find(sdiff==1);
soff = find(sdiff==-1);

cs = colormap('cool');
cs = cs(round(linspace(36,64,5)),:);
cs = max(0, cs-.1);
rng(1);
cs = cs(randperm(size(cs,1)),:);

yh=.02;
X = [.05 .2 .35 .72];
Y = [.05 .05+dy*.55 .64];

clf;

clear hs;
i=0;

i=i+1;
hs{i} = my_subplot(2,3,1,[.95 .85]);
hs{i}.Position(1) = hs{i}.Position(1) +.04;
hs{i}.Position(2) = hs{i}.Position(2) +.01;
hold all;
cpc = colormap('winter');
cpc = cpc(10:6:end,:);
for j = 1:10
    fj = faceEx(:,j);
    fj = fj * sign(skewness(fj));
    fj = fj(tpts);
    fj = bin2d(fj(:),4);
    plot(tpts(1:4:end),my_conv2(fj,4,1)/nanstd(fj)+(10-j)*2.5+2,'color',cpc(j,:),'linewidth',.5);
end
ht=text(-.08,.75,'Face motion energy PCs','fontsize',8,'fontangle','normal','color',cpc(5,:),...
    'HorizontalAlignment','center');
ht.Rotation=90;
ineur = 1:10:10*10;
for j = 1:10 %size(Fex,2)
    fplot = Fex(ineur(j),tpts);
    fplot = bin2d(fplot(:),4);
    plot(tpts(1:4:end),zscore(fplot)/4-j*2.5,'color','k','linewidth',.5)%cn(j,:))
end
axis off;
axis tight;
plot((find(diff(stimtpt(tpts))==1,1)+tpts(1)-1) * [1 1], [-10*2.5 10.5*2.5],'k--','color',.5*[1 1 1]);
plot((find(diff(stimtpt(tpts))==-1,1)+tpts(1)-1) * [1 1], [-10*2.5 10.5*2.5],'k--','color',.5*[1 1 1]);
ht=text(-.08,.25,'Example neurons','fontsize',8,'fontangle','normal',...
    'HorizontalAlignment','center');
ht.Rotation=90;

axes('position',[hs{i}.Position(1) hs{i}.Position(2)-yh hs{i}.Position(3) yh*.8]);
fig_stimon(stimtpt(tpts))

cm = colormap('hot');
cm = cm(32+[1:4:4*4],:);

i=i+1;
hs{i} = my_subplot(6,5,3,[.55 .7]);
loglog([1e-3 1],[1e-3 1],'k','linewidth',.5)
hold all;
for j = 1:4
    x = vface(1:10,1,j)/nansum(vface(:,1,j))*1;
    y = vface(1:10,2,j)/nansum(vface(:,2,j))*1;
    loglog(x,y,'.','color',cm(j,:),'markersize',8)
end
%ylim([.8 1.2]);
box off;
axis tight;
set(gca,'ytick',10.^[-2:0],'yticklabel',{'0.01','0.1','1'},...
    'xtick',10.^[-2:0],'xticklabel',{'0.01','0.1','1'});
ylabel({'spont periods'});
xlabel('stim periods');
%title('% variance of PCs','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

i=i+1;
hs{i} = my_subplot(6,5,4,[.55 .7]);
hs{i}.Position(1) = hs{i}.Position(1) - .0;
%hs{i}.Position(1) = hs{i-1}.Position(1) + hs{i-1}.Position(3)*1.3;
hold all;
bed = [-.6:.1:10];
for d = 1:4
    histogram(snr{d},'binedges',bed,'displaystyle','stairs',...
        'normalization','probability','edgecolor',cm(d,:));
    plot(nanmean(snr{d}),.07,'v','color',cm(d,:),...
        'markersize',4);
end
axis tight;
xlabel('tuning SNR');
ylabel('fraction');
axis square;


dex = 1;
[sval,sbins] = sort(sv{dex});
[fval,fbins] = sort(facevar{dex});

i = i+1;
hs{i}=my_subplot(6,5,5,[.55 .7]);
hs{i}.Position(1) = hs{i}.Position(1) - .0;
plot(sbins(1:3:end),fbins(1:3:end),'.','markersize',2);
xlabel({'from stimulus (rank)'});
ylabel({'from behavior (rank)'});
box off;
axis tight;
axis([1 max(sbins) 1 max(sbins)]);
axis square;
set(gca,'ytick',10.^[0 4],'yticklabel',{'1','10^4'},...
    'xtick',10.^[0 4],'xticklabel',{'1','10^4'});

i = i+1;
hs{i}=my_subplot(6,5,8,[.55 .7]);
axis square;
axis off;
hp=hs{i}.Position;
hp(2)=hp(2)-.085;
hp(4)=hp(4)+.06;
hp(3)=.07;
hss=axes('position',hp);
ss = sface{dex}(iface{dex},:);
ss = my_conv2(ss,10,1);
ss = my_conv2(ss,2,1);
imagesc(ss,[-.5 .5])
set(gca,'ytick',10.^[0 4],'yticklabel',{'1','10^4'},...
    'xtick',[]);
box off;
text(.5,-.05,{'behavior','coefficients'},'HorizontalAlignment','center','fontsize',8,'fontangle','normal');
ylabel('behavior embedding rank');
colormap(hss,redblue);

hp2=hp;
hp2(1)=hp(1)+hp(3)+.025;
hp2(3)=hp2(4);
axes('position',hp2);
plot(istim{dex}(1:3:end),iface{dex}(1:3:end),'.','markersize',2);
box off;
axis tight;
axis([1 max(sbins) 1 max(sbins)]);
axis square;
set(gca,'ytick',10.^[0 4],'yticklabel',{'1','10^4'},...
    'xtick',10.^[0 4],'xticklabel',{'1','10^4'},...
    'ydir','reverse');

hp3=hp2;
hp3(4)=hp(3);
hp3(2)=hp3(2)-hp3(4)-.025;
axes('position',hp3);
ss = sstim{dex}(istim{dex},:);
ss = my_conv2(ss,10,1);
[~,imax] = max(ss,[],1);
[~,isort] = sort(imax);
ss = sstim{dex}(istim{dex},:);
ss = my_conv2(ss,3,1);
imagesc(ss(:,isort(end:-1:1))',[-.5 .5])
set(gca,'xtick',10.^[0 4],'xticklabel',{'1','10^4'},...
    'ytick',[]);
box off;
ht=text(-.11,.5,'stimuli','HorizontalAlignment','center','fontsize',8,'fontangle','normal');
ht.Rotation=90;
xlabel('stimulus embedding rank');
colormap(gca,redblue);


i = i+1;
hs{i}=my_subplot(6,5,10,[.55 .7]);
axis off;
hp=axes('position',hp2);
hp.Position(1)=hp.Position(1)+.25;
hs{i}.Position(1)=hp.Position(1);
pds = abs(istim{dex}(1:7:end) - istim{dex}(1:7:end)');
pds = pds - diag(NaN*diag(pds));
pdf = abs(iface{dex}(1:7:end) - iface{dex}(1:7:end)');
pdf = pdf - diag(NaN*diag(pdf));
plot(pds(1:501:end), pdf(1:501:end), '.','markersize',2);
set(gca,'ytick',10.^[0 4],'yticklabel',{'1','10^4'},...
    'xtick',10.^[0 4],'xticklabel',{'1','10^4'});
p = polyfit(pds(~isnan(pds)),pdf(~isnan(pds)),1);
hold all;
plot([0 1.2e4],[0 1.2e4] * p(1) + p(2),'k');
[r,p]=corr(pds(~isnan(pds)),pdf(~isnan(pds)));
text(.65,.45,sprintf('r = %1.3f',r),'fontsize',8,'fontweight','bold','fontangle','normal');
ylabel({'from behavior'})
xlabel({'from stimulus'})
axis square;
axis tight;
box off;


i=i+1;
hs{i} = my_subplot(6,5,16,[.55 .7]);
hs{i}.Position(1) = hs{i}.Position(1) +.01;
hs{i}.Position(2) = hs{i}.Position(2) +.0;
hold all;
for j = 1:size(Vshared,2)
    plot(Vshared(:,j,1)*100,'.-','markersize',12,'color',cm(j,:),'linewidth',.5)
end
box off;
ylabel({'% stim variance'});
xlabel('dimension');
%title('stim-face subspace','fontweight','normal','fontangle','italic','fontsize',10)
axis square;


i=i+1;
hs{i} = my_subplot(6,5,17,[.55 .7]);
hs{i}.Position(2) = hs{i}.Position(2) +.0;
hold all;
for j = 1:size(Ushared,1)
    uplot = Ushared{j,kt};
    uplot = uplot * sign(mean(uplot));
    bed = [-.02:.001:.035];
    histogram(uplot,20,'binedges',bed,'edgecolor',cm(j,:),...
        'normalization','probability','displaystyle','stairs');
end
box off;
xlim([-.02 .035]);
ylabel('fraction');
xlabel('neural weights');
%title('1st PC weights','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

i=i+1;
hs{i} = my_subplot(3,2,5,[.7 .5]);
hs{i}.Position(1)=hs{i-2}.Position(1);
axis off;
hp=my_subplot(3,2,5,[.9 .9]);
hp.Position(1)=hp.Position(1)-.05;
hp.Position(2)=hp.Position(2)-.01;
im=imread('schematicNEW.png');
image(im);
axis tight;
axis image
axis off;

kt = 1; % use face projections
tshared1 = tshared{kt} * sign(mean(Ushared{1,kt}));
tshared2 = tshared{2} * sign(mean(Ushared{2,2}));
i=i+1;
hs{i} = my_subplot(2,2,4,[.9 .65]);
hs{i}.Position(1) = hs{5}.Position(1);
hs{i}.Position(2) = hs{i}.Position(2)-.01;
hs{i}.Position(3) = .5;
dobin=1;
tproj1=tprojS{1}(1:3,:);
tproj1(1,:) = tproj1(1,:)*-1;
tproj1(3,:) = tproj1(3,:)*-1;
tproj2=tprojS{2}(1:3,:);
fig_traces(sprojS{kt}(1:5,tpts),tshared1(tpts),tproj1(:,tpts),tshared2(tpts),tproj2(:,tpts),cs,dobin);
axes('position',[hs{i}.Position(1) hs{i}.Position(2)-yh hs{i}.Position(3) yh*.8]);
fig_stimon(stimtpt(tpts),0);
axes('position',[hs{i}.Position(1) hs{i}.Position(2)-4*yh hs{i}.Position(3) yh*2.8]);
plot(bin2d(runspeed(tpts),4),'color',[.2 .8 .2]);
text(.2,1.0,'running','color',[.2 .8 .2],'fontangle','normal','fontsize',8);
hold all;
plot(1 + [0 60*5],-6 * [1 1],'k','linewidth',2);
text(0.06,-.1,'5 min','horizontalalignment','center','fontsize',8,'fontweight','bold','fontangle','normal');
axis off;
axis tight;


tstr{1} = 'Behavior and neural activity with/without visual stimulation';
tstr{2} = 'Variance of face PCs';
tstr{3} = 'Neural stimulus tuning';
tstr{4} = 'Variance explained';
tstr{5} = 'Comparison of embeddings';
tstr{6} = 'Embedding distances';
tstr{7} = 'Stim-face dimensions';
tstr{8} = 'Top dimension';
tstr{9} = '';
tstr{10} = 'Projections of neural activity';

% -------------- LETTERS
hp=.05;
hy=1.22;
deffont=8;
for j = [1:length(hs)]
    if j ==1
        hp0 =.04;
        hy0 = 1.06;
    elseif j==5 || j==6 || j==length(hs)
        hy0=1.08;
    else
        hp0=hp;
        hy0 = hy;
    end
    hpos = hs{j}.Position;
    lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
    axes('position', lpos);
    text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
    
    axes('position', [lpos(1)+.02 lpos(2)-.005 lpos(3) lpos(4)]);
    text(0,0, tstr{j},'fontangle','normal','fontsize',deffont);
    axis([0 1 0 1]);
    axis off;
end
