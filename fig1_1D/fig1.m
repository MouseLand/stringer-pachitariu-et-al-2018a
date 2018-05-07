function fig1(matroot)

load(fullfile(matroot,'corr1stpc.mat'));

%%
default_figure([1 1 6.25 8]);

%%
%trange = 750+[1:1500];
trange = 3000+[1:2000];

spkplot=(my_conv2(zscore(results.spks(end:-1:1,trange),1,2),[10 1],[1 2]));
%%
clear hs;
i = 0;
cm = colormap('gray');

clf;

ndat = numel(results.runcorr);
load('cdat.mat');
xh = .5;
yh = .5;
yd1=.12;
yd2 = .07;

i=i+1;
hs{i}=my_subplot(5,4,1,[xh yh]);
axis off;
hp=hs{i}.Position;
axes('position',[hp(1)-.05 hp(2)-.03 1.4*hp(3:4)]);
im = imread('planes.PNG');
imagesc(im);
axis off;
axis image;

i=i+1;
hs{i}=my_subplot(5,4,2,[xh yh]);
axis off;
hp=hs{i}.Position;
axes('position',[hp(1)-.05 hp(2)-.03 1.4*hp(3:4)]);
im = imread('explane.PNG');
imagesc(im);
title('plane 5/11','fontweight','normal','fontsize',8);
axis off;
axis image;

% ---------- correlation histogram -------------- %
i = i+1;
hs{i} = my_subplot(5,4,3,[xh yh]);
bins = [-.3:.001:.75];
nb = histcounts(results.cshuff(1:1:end),bins);
%nb = nb/sum(nb);
%area(bins(1:end-1),log(nb),-17,'facecolor',.5*[1 1 1],'edgecolor','none','showbaseline','off')
pbins = bins(1:end-1)+(bins(2)-bins(1))/2;
plot(pbins, log10(nb), 'color', .5*[1 1 1]);
hold all;
nb = histcounts(results.ccall(1:1:end),bins);
%nb = nb/sum(nb);
plot(pbins,log10(nb),'color',cdat(1,:))
text(.5,.75,'shuffled','fontangle','normal','fontsize',8,'color',.5*[1 1 1],'fontweight','bold');
text(.5,.88,'original','fontangle','normal','fontsize',8,'color',cdat(1,:),'fontweight','bold');
axis tight;
%ylim([-16 -3]);
xlim([-.25 .7]);
set(gca,'ytick',[0:2:6],'yticklabel',{'10^0','10^2','10^4',...
    '10^6'});
xlabel('correlation');
ylabel('number of pairs');
box off;
axis square;

% ---------- mean and std of correlations -------------- %
i = i+1;
hs{i} = my_subplot(5,4,4,[xh yh]);
hold all;
id = [1 2 5 3 4 6:9]; % show gray screen first
for d = 1:ndat
    errorbar(d, results.mcorr(id(d)), results.stdcorr(id(d),1), results.stdcorr(id(d),2),'.', 'color', cdat(id(d),:))
end
hi = .23;
y1 = [.5 3.7];
y2 = [3.3 9.5];
for j = 1:2
plot([y1(j) y2(j)],hi*[1 1],'k');
plot(y1(j)*[1 1],[hi-.02 hi],'k');
plot(y2(j)*[1 1],[hi-.02 hi],'k');
end
text(0.02,1.26,{' gray','screen'},'fontsize',7,'fontangle','normal')
text(0.45,1.11,{'darkness'},'fontsize',7,'fontangle','normal')
ylabel('correlation');
xlabel('recordings');
set(gca,'xtick',[]);
box off;
axis square;
ylim([-.05 .25]);

%
% -------------------- CORR MATRICES ---------%
ineu1 = 1:20:20*10;
ineu2 = 2:20:20*10;
rng(6);
ineu1 = randperm(size(results.cc1,1),10);
ineu1 = randperm(size(results.cc1,1),10);
crange = [-.15 .15];
i = i+1;
hs{i} = my_subplot(5,4,5,[xh yh]);
axis off;
hp=hs{i}.Position;
hp(1)=hp(1)-.03;
hp(2)=hp(2)-.01;
hp(3)=hp(3)+.025;
hp(4)=hp(4)+.02;
axes('position',hp);
imagesc(results.cc1(ineu1,ineu2), crange);
text(.5, 0,'neurons 11-20','HorizontalAlignment','center','fontangle','normal','fontsize',8);
ht=text(-.17, .5,'neurons 1-10','HorizontalAlignment','center','fontangle','normal','fontsize',8);
ht.Rotation=90;
axis off;
axis square;
colormap(gca,'redblue');
axes('position',[hp(1)+hp(3) hp(2) .03 hp(4)]);
cb=colorbar;
cb.Ticks=[1/6 .5 5/6];
cb.TickLabels={'-0.1','0','0.1'};
colormap(gca,'redblue');
axis off;

i = i+1;
hs{i} = my_subplot(5,4,6,[xh yh]);
axis off;
hp=hs{i}.Position;
hp(1)=hp(1)-.03;
hp(2)=hp(2)-.01;
hp(3)=hp(3)+.025;
hp(4)=hp(4)+.02;
axes('position',hp);
imagesc(results.cc2(ineu1,ineu2), crange);
axis off;
axis square;
colormap(gca,'redblue');
text(.5, 0,'neurons 11-20','HorizontalAlignment','center','fontangle','normal','fontsize',8);
ht=text(-.17, .5,'neurons 1-10','HorizontalAlignment','center','fontangle','normal','fontsize',8);
ht.Rotation=90;

% --- corr halves X Y ---- %
i = i+1;
hs{i} = my_subplot(5,4,7,[xh yh]);
plot(results.cc1(1:50:end), results.cc2(1:50:end), '.','markersize',4,...
    'color',cdat(1,:));
axis tight;
box off;
xlabel('correlation half 1');
ylabel('correlation half 2');
text(.25,1.05, sprintf('r=%1.2f',results.rcorr(dex)));
axis square;

% ------------ R OF CORRELATIONS
i = i+1;
hs{i} = my_subplot(5,4,8,[xh yh]);
histogram(results.rcorr,'binedges',[.5:.1:.9],'edgecolor',[0 0 0],...
    'facecolor', .1*[1 1 1]);
box off;
ylabel('number of recordings');
xlabel('r of correlations');
xlim([0 1]);
axis square;


% --------- SPIKE plot ------------------ %
i = i+1;
x0 = .26;
y0 = .76;
xh0 = .7;
yh0 = .24;
hs{i} = my_subplot(4,1,3,[.86 .78]);%axes('position',[x0 y0 xh0 yh0]);
hs{i}.Position(2)=hs{i}.Position(2)+.12;
hs{i}.Position(1)=hs{i-4}.Position(1);
pos=hs{i}.Position;
[NN,NT] = size(results.spks);
imagesc(spkplot(:,:), [0 .3]);
hold all;
colormap(hs{i},flipud(cm));
plot([1 1]+10+size(spkplot,2), NN-1000+[0 1000],'k','linewidth',2);
ht=text(1.03,-.0,'1000 neurons','fontangle','normal','fontsize',8,'HorizontalAlignment','right');
ht.Rotation = 270;
%ht=text(-0.03,.5,'sorted by 1st PC weights','fontangle','normal','fontsize',8,...
%    'HorizontalAlignment','center','fontweight','bold');
%ht.Rotation = 90;
axis off;
axis tight;

axes('position',[pos(1) pos(2)-.15 pos(3) .14]);
hold all;
pcs = results.firstpcs(trange,:);
pcs = pcs - mean(pcs,1);
pcs = pcs / std(pcs(:,1));
j=1;
cm(1,:) = [.2 .8 .2];
beh = min(5,zscore(results.behavior(trange(2:numel(trange)),j)));
plot(beh,'linewidth',1,'color',cm(j,:),'linewidth',1);
text(.0,1.02,'running','fontangle','normal','fontsize',8,'color',[.2 .8 .2],...
    'fontweight','bold');

pstr={'1st PC','2nd PC','3rd PC'};
spc = sign(skewness(pcs(:,1:3)));
for j = 1:3
    plot(pcs(:,j)*spc(j) - (j-1)*5,'k','linewidth',.5);
    text(-5,-(j-1)*5+1, pstr{j},'units','data','fontangle','normal','fontsize',8,...
        'HorizontalAlignment','right');
end
plot([0 300/1.2],[-1 -1]*(j)*5+3,'k')
text(0,0,'  5 minutes','fontangle','normal','fontsize',8);
axis tight;
axis off;
%


% ------ 1st PC stats ---------------- %
% 
% i = i+1;
% hs{i} = my_subplot(3,5,11,[xh yh]);
% hs{i}.Position(2) = hs{i}.Position(2) -yd2;
% histogram(results.pcdist,'binedges',[.2:.1:.8],'edgecolor',[0 0 0],...
%     'facecolor', .1*[1 1 1]);
% box off;
% ylabel('number of recordings');
% xlabel('% of positive weights');
% xlim([0 1]);
% axis square;

% ----------- correlation with behavior ----------%
% cm is running pupil whisker
cm(1,:) = [.2 .8 .2];
cm(2,:) = [.5 .6 .5];
cm(3,:) = [0 .2 0];
tstr = {'running','pupil area','whisking'};
i = i+1;
hs{i} = my_subplot(5,4,17,[xh yh]);
hold all; 
cbeh = abs([results.runcorr(:) results.pupilcorr(:) results.whiskcorr(:)]);
for j = 1:3
    for d = 1:ndat
        plot(j+randn*.1,cbeh(d,j),'.','color',cdat(d,:),'markersize',6);
    end
    errorbar(j,mean(cbeh(:,j)),std(cbeh(:,j))/sqrt(ndat-1),'k.','linewidth',1,'markersize',8);
    
end
%plot(abs(results.runcorr),abs(results.pupilcorr),'o','color',cm(2,:),'markersize',4);
axis square;
%axis tight;
axis([.5 3.5 0 1]);
set(gca,'xtick',[1 2 3],'xticklabel',{'running','pupil area','whisking'});
xtickangle(45);
ylabel('correlation with 1st PC');


% --------- correlation vs distance
i = i+1;
hs{i} = my_subplot(5,4,18,[xh yh]);
dx = results.dbins(1:end-1) + diff(results.dbins(1:2));
hold all;
for j = [1:size(results.dmean,2)]
    plot(dx, results.dmean(:,j),'color',cdat(j,:),'linewidth',.25)
end
j=dex;
plot(dx, results.dmean(:,j),'color',cdat(j,:),'linewidth',1)
xlabel('distance (um)');
ylabel('mean correlation');
ylim([-.1 .1]);
axis square;

% ----------- example plane -------%
i = i+1;
hs{i} = my_subplot(5,4,19,[xh yh]);
axis off;
hp=hs{i}.Position;
hp(1)=hp(1)-.045;
hp(2)=hp(2)-.015;
hp(3)=hp(3)+.035;
hp(4)=hp(4)+.035;
axes('position',hp);
iplane = 35*5;
icell = find(results.cellmed(:,3)==iplane);
med = results.cellmed(icell,1:2);
ctype= results.cellpc(icell);
hold all;
cm=[1 0 0; 0 0 1];
NN = size(med,1);
for j = randperm(NN)
    plot(med(j,1),med(j,2),'.','markersize',3,'color',cm(ctype(j),:));
end
axis off;
axis square;
plot([0 500],[0 0]-30,'k','linewidth',2)
text(0,0,'500 \mum','fontsize',8,'fontangle','normal');
text(.45, 0, {'positive weight'},'color',cm(1,:),'fontsize',8,'fontangle','normal');
text(.45, -.1, {'negative weight'},'color',cm(2,:),'fontsize',8,'fontangle','normal');

i = i+1;
hs{i} = my_subplot(5,4,20,[xh yh]);
hold all;
for d = 1:ndat
    lw = .5;
    if d==dex
        lw = 1;
    end
    plot([results.indist(:,d) results.outdist(:,d)]','color',cdat(d,:),'linewidth',lw)
end
p = signrank(results.indist(:),results.outdist(:));
ylim([0 580]);
xlim([.75 2.25]);
box off;
ylabel('pairwise distance (um)');
set(gca,'xtick',[1 2],'xticklabel',{'same sign',['opp sign']});
axis square;

% %
% for j = [1 2 3]
%     switch j
%         case 1
%             cc = abs(results.runcorr);
%         case 2
%             cc = abs(results.pupilcorr);
%         case 3
%             cc = abs(results.whiskcorr);
%     end
%     
%     histogram(cc,'binedges',[0:.05:1],'facecolor','none',...
%         'edgecolor',cm(j,:),'displaystyle','stairs','linewidth',2);
%     box off;
%     ylabel('number of recordings');
%     xlabel('correlation with 1st PC');
%     %title('1st PC correlation with running')
%     text(.03,1-.13*(j-1),tstr{j},'color',cm(j,:),'fontangle','normal','fontsize',8);
%     xlim([0 1]);
%     axis square;
%     axis tight;
% end

%
tstr{1} = 'imaging setup';
tstr{2} = '';
tstr{3} = '';
tstr{4} = '';
tstr{5} = 'correlations half 1';
tstr{6} = 'correlations half 2';
tstr{7} = '';
tstr{8} = {''};%correlation of halves'};
tstr{11} ='';
tstr{10} ={''};
tstr{9} ='Population activity sorted by 1st PC weights';
tstr{12} =  '';
tstr{13} = '';%'Spatial distributions';

% -------------- LETTERS
hp=.09;
hy=1.4;
deffont=8;
for j = [1:length(hs)]
    if j==9
        hp0 =hp;
        hy0 = 1.18;
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
        
    if j==10
        tsh=.02;
    else
        tsh=.03;
    end
    axes('position', [lpos(1)+tsh lpos(2)-.006 lpos(3) lpos(4)]);
    text(0,0, tstr{j},'fontangle','normal','fontsize',8);
    axis([0 1 0 1]);
    axis off;
    
end



