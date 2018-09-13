function fig1new(matroot)

pc = load(fullfile(matroot,'corr1stpc.mat'));
clust = load(fullfile(matroot,'clust1D.mat'));
try
	dd = load(fullfile(matroot,'PCApred.mat'));
catch
    dd = load('PCApred.mat');
end
expv_neurons = dd.expv_neurons;
ndims0 = 2.^[0:size(expv_neurons,1)-1];

dbeh=load(fullfile(matroot,'expv_behavior_neurons.mat'));
expv_behavior = dbeh.expv_behavior;

%%
default_figure([1 1 6.25 7]);

%%
%trange = 750+[1:1500];
trange = [1:400];

pcplot=(my_conv2(zscore(pc.results.spks(end:-2:1,trange),1,2),[3 1],[1 2]));
clustplot=(my_conv2(zscore(clust.results.spks(end:-1:1,trange(1):trange(end)*3),1,2),[5 1],[1 2]));
%%
clear hs;
i = 0;
cm = colormap('gray');


% PC colors
nPCs = 4;
cpc = max(0,colormap('spring'));
cpc = cpc(1:32,:);
cpc = cpc(round(linspace(1,size(cpc,1),nPCs)),:);

clf;

ndat = numel(pc.results.runcorr);
load('cdat.mat');
xh = .55;
yh = .55;
yd1=.12;
yd2 = .07;

i=i+1;
hs{i}=my_subplot(5,4,1,[xh yh]);
axis off;
hp=hs{i}.Position;
axes('position',[hp(1)-.05 hp(2)-.05 1.4*hp(3:4)]);
im = imread('planes.PNG');
imagesc(im);
axis off;
axis image;

i=i+1;
hs{i}=my_subplot(5,4,2,[xh yh]);
axis off;
hp=hs{i}.Position;
axes('position',[hp(1)-.05 hp(2)-.05 1.4*hp(3:4)]);
im = imread('explane.PNG');
imagesc(im);
title('plane 5/11','fontweight','normal','fontsize',8);
axis off;
axis image;

% ---------- correlation histogram -------------- %
i = i+1;
hs{i} = my_subplot(5,4,3,[xh yh]);
bins = [-.3:.001:.75];
nb = histcounts(pc.results.cshuff(1:1:end),bins);
%nb = nb/sum(nb);
%area(bins(1:end-1),log(nb),-17,'facecolor',.5*[1 1 1],'edgecolor','none','showbaseline','off')
pbins = bins(1:end-1)+(bins(2)-bins(1))/2;
plot(pbins, log10(nb), 'color', .5*[1 1 1]);
hold all;
nb = histcounts(pc.results.ccall(1:1:end),bins);
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
    errorbar(d, pc.results.mcorr(id(d)), pc.results.stdcorr(id(d),1), pc.results.stdcorr(id(d),2),'.', 'color', cdat(id(d),:))
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


% --------- PC img ------------------ %
i = i+1;
x0 = .26;
y0 = .76;
xh0 = .7;
yh0 = .24;
hs{i} = my_subplot(4,2,3,[1.3 .9]);%axes('position',[x0 y0 xh0 yh0]);
hs{i}.Position(2)=hs{i}.Position(2)+0.02;
hs{i}.Position(1)=hs{i}.Position(1)+.14;
%hs{i}.Position(1)=hs{i-4}.Position(1);
pos=hs{i}.Position;
[NN,NT] = size(pc.results.spks);
imagesc(pcplot(:,:), [0 .3]);
hold all;
colormap(hs{i},flipud(cm));
%plot([1 1]+10+size(pcplot,2), NN-1000+[0 1000],'k','linewidth',2);
ht=text(1.03,-.0,'1000 neurons','fontangle','normal','fontsize',8,'HorizontalAlignment','right');
ht.Rotation = 270;
%ht=text(-0.03,.5,'sorted by 1st PC weights','fontangle','normal','fontsize',8,...
%    'HorizontalAlignment','center','fontweight','bold');
%ht.Rotation = 90;
axis off;
axis tight;

% ---------- CLUST img ------------- %
i=i+1;
hs{i} = my_subplot(4,2,5,[1.1 0.9]);%axes('position',[x0 y0 xh0 yh0]);
%hs{i}.Position(2)=hs{i}.Position(2)+.1;
hs{i}.Position(1) = hs{i-1}.Position(1);
hs{i}.Position(4) = hs{i-1}.Position(4);
hs{i}.Position(3) = hs{i-1}.Position(3);
pos=hs{i}.Position;
[NN,NT] = size(clustplot);
imagesc(clustplot(:,:), [0 .3]);
hold all;
colormap(hs{i},flipud(cm));
%plot([1 1]+10+size(clustplot,2), NN-1000+[0 1000],'k','linewidth',2);
ht=text(1.03,-.0,'1000 neurons','fontangle','normal','fontsize',8,'HorizontalAlignment','right');
ht.Rotation = 270;
%ht=text(-0.03,.5,'sorted by 1st PC weights','fontangle','normal','fontsize',8,...
%    'HorizontalAlignment','center','fontweight','bold');
%ht.Rotation = 90;
axis off;
axis tight;

axes('position',[pos(1) pos(2)-.13 pos(3) .1]);
hold all;
cm(1,:) = [.2 .8 .2];
cm(2,:) = [.5 .6 .5];
cm(3,:) = [0 .2 0];
tstr = {'running','pupil area','whisking'};
for j = 1:3
	beh = min(5,zscore(pc.results.behavior(trange(2:numel(trange)),j)));
	if j==3
		beh=beh*-1;
	end
	plot(beh,'linewidth',1,'color',cm(j,:),'linewidth',1);
	text(.0,1.15-(j-1)*0.18,tstr{j},'fontangle','normal','fontsize',8,'color',cm(j,:),...
		'fontweight','bold');
end
axis tight;
axis off;

axes('position',[pos(1) pos(2)-.23 pos(3) .1]);
pcs = pc.results.firstpcs(trange,:);
pcs = pcs - nanmean(pcs,1);
pcs = pcs / nanstd(pcs(:,1));
pstr={'PC1','PC2','PC3'};
spc = sign(skewness(pcs(:,1:3)));
for j = 1:3
    plot(pcs(:,j)*spc(j) - (j-1)*3,'k','linewidth',.5);
    text(-5,-(j-1)*3+1.2, pstr{j},'units','data','fontangle','normal','fontsize',8,...
        'HorizontalAlignment','right');
	hold all;
end
%plot([0 300/1.2],[-1 -1]*(j)*5+3,'k')
text(0,0,'  5 minutes','fontangle','normal','fontsize',8);
axis tight;
axis off;

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
tstr = {'running','pupil area','whisking'};
i = i+1;
hs{i} = my_subplot(5,4,8,[xh yh]);
%hs{i}.Position(2) = hs{i}.Position(2)+0.05;
hold all; 
cbeh = abs([pc.results.runcorr(:) pc.results.pupilcorr(:) pc.results.whiskcorr(:)]);
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

% -------- VAR EXP BY BEH --------------- %
i=i+1;
hs{i} = my_subplot(5,4,12,[xh yh]);
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

% ------ PEER PRED DIAGRAM ------ %
i=i+1;
hs{i}=my_subplot(5,4,16,[xh yh]);
axis off;
hold all;
msize = 6;
nls = [8 4 1];
for j = 1:length(nls)
    call{j} = repmat([0 0 0],nls(j),1);
end
call{2} = cpc;
ells = [3 3 0];
hp=hs{i}.Position;
hp(1)=hp(1)-.05;
hp(3)=hp(3)+.03;
hp(2)=hp(2)-.04;
hp(4)=hp(4)+.02;
layeredNet(hp, nls, call, msize,ells);
text(0,1.23,'peers','HorizontalAlignment','center','fontsize',8);
text(.5,1.03,'neural PCs','HorizontalAlignment','center','fontsize',8);
text(1,0.5,{'single','neuron','activity'},'HorizontalAlignment','center','fontsize',8);
axis off;

% ----- PEER PRED ---------%
i=i+1;
hs{i}=my_subplot(5,4,20,[xh yh]);
for d = 1:ndat
    semilogx(ndims0,expv_neurons(:,d),'color',cdat(d,:),'linewidth',1);
    hold all;
end
semilogx(ndims0,mean(expv_neurons,2),'k','linewidth',2);
box off;
axis square;
xlabel('dimensions');
ylabel({'variance explained','(test set)'});
grid on;
grid minor;
grid minor;
set(gca,'xtick',2.^[0:4:10],'ytick',[0:.10:.30]);



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
tstr{5} = '';
tstr{6} = '';
tstr{7} = '';
tstr{8} = {''};%correlation of halves'};
tstr{11} ='';
tstr{10} ={''};
tstr{9} ='';
tstr{12} =  '';
tstr{13} = '';%'Spatial distributions';

% -------------- LETTERS
hp=.09;
hy=1.4;
deffont=8;
for j = [1:length(hs)]
    if j==5
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



