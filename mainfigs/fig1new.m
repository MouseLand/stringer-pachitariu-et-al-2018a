function fig1new(matroot)

pc = load(fullfile(matroot,'corr1stpc.mat'));
clust = load(fullfile(matroot,'clust1D.mat'));

dbeh=load(fullfile(matroot,'expv_behavior_neurons.mat'));
expvPC_behavior = dbeh.expvPC_behavior;

load(fullfile(matroot,'PCpredict.mat'));

%%
default_figure([1 1 6.25 7.5]);

%%
%trange = 750+[1:1500];
trange = [1:600];

pcplot=(my_conv2(zscore(pc.results.spks(1:2:end,trange),1,2),[3 0.5],[1 2]));
clustplot=(my_conv2(zscore(clust.results.spks(end:-1:1,trange),1,2),[3 0.5],[1 2]));
%%
clear hs;
i = 0;
cgray = colormap('gray');


% PC colors
nPCs = 20;
cpc = max(0,colormap('spring'));
cpc = cpc(1:48,:);
cpc = cpc(round(linspace(1,size(cpc,1),nPCs)),:);

clf;

ndat = numel(pc.results.runcorr);
load('cdat.mat');
xh = .5;
yh = .5;
yd1=.12;
yd2 = .07;

i=i+1;
hs{i}=my_subplot(5,5,1,[xh yh]);
axis off;
hp=hs{i}.Position;
axes('position',[hp(1)-.05 hp(2)-.05 1.4*hp(3:4)]);
im = imread('planes.PNG');
imagesc(im);
axis off;
axis image;

i=i+1;
hs{i}=my_subplot(5,5,2,[xh yh]);
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
hs{i} = my_subplot(5,5,3,[xh yh]);
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
hs{i} = my_subplot(5,5,4,[xh yh]);
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

pcs = pc.results.firstpcs(trange(1:length(trange)/3),:);
pcs = zscore(pcs.*sign(skewness(pcs)),1,1);
pcs = pcs - min(pcs,[],1);
behall = min(5,zscore(pc.results.behavior(trange(1:length(trange)/3),:)));
behall(:,3) = behall(:,3)*-1;
behall = behall - min(behall,[],1);

% ----------- correlation with behavior ----------%
tstr = {'running','pupil area','whisking','all 3'};
i = i+1;
hs{i} = my_subplot(5,5,5,[xh yh]);
%hs{i}.Position(2) = hs{i}.Position(2)+0.05;
hold all; 
plot(pcs(:,1),behall(:,1),'.','color',cpc(j,:));
xlabel('PC magnitude');
ylabel({'running','speed (a.u.)'});
axis square;
axis tight;
title({'PC 1','r=0.75'},'fontweight','normal');

% --------- PC img ------------------ %
i = i+1;
x0 = .26;
y0 = .76;
xh0 = .7;
yh0 = .24;
hs{i} = my_subplot(5,1,2,[.92 .8]);%axes('position',[x0 y0 xh0 yh0]);
hs{i}.Position(2)=hs{i}.Position(2)-0.1;
hs{i}.Position(1)=hs{i}.Position(1);
%hs{i}.Position(1)=hs{i-4}.Position(1);
pos=hs{i}.Position;
[NN,NT] = size(pcplot);
imagesc(pcplot(:,:), [0 .3]);
hold all;
colormap(hs{i},flipud(cgray));
plot(-5*[1 1], [NN-500 NN],'k','linewidth',2);
ht=text(-.03,-.0,'1000 neurons','fontangle','normal','fontsize',8,'HorizontalAlignment','left');
ht.Rotation = 90;
%ht=text(-0.03,.5,'sorted by 1st PC weights','fontangle','normal','fontsize',8,...
%    'HorizontalAlignment','center','fontweight','bold');
%ht.Rotation = 90;
axis off;
axis tight;

axes('position',[pos(1) pos(2)+pos(4)+.02 pos(3) .08]);
hold all;
cm(1,:) = [.2 .8 .2];
cm(2,:) = [.5 .6 .5];
cm(3,:) = [.1 .9 .6];
tstr = {'running','pupil area','whisking'};
for j = 1:3
	beh = behall(:,j);
	plot(beh,'linewidth',1,'color',cm(j,:),'linewidth',1);
	text(.0,1.1-(j-1)*0.18,tstr{j},'fontangle','normal','fontsize',8,'color',cm(j,:));
		%'fontweight','bold');
end
plot(pcs(:,1),'-.','color',cpc(1,:),'linewidth',0.5);
text(.2,1,'PC 1','HorizontalAlignment','right','fontangle','normal','color',cpc(1,:),'fontsize',8)
text(0,-.05,'  10 s','fontangle','normal','fontsize',8);
plot([0 10/1.2],-.3*[1 1],'k','LineWidth',2);
axis tight;
axis off;
xlim([-5/3 size(pcs,1)]);

% ---------- CLUST img ------------- %
i=i+1;
hs{i} = my_subplot(5,1,3,[1 0.9]);%axes('position',[x0 y0 xh0 yh0]);
hs{i}.Position(2)=hs{i}.Position(2)-.06;
hs{i}.Position(1) = hs{i-1}.Position(1);
hs{i}.Position(4) = hs{i-1}.Position(4);
hs{i}.Position(3) = hs{i-1}.Position(3);
pos=hs{i}.Position;
[NN,NT] = size(clustplot);
imagesc(clustplot(:,:), [0 .3]);
hold all;
colormap(hs{i},flipud(cgray));
%plot([1 1]+10+size(clustplot,2), NN-1000+[0 1000],'k','linewidth',2);
plot(-5*[1 1], [NN-500 NN],'k','linewidth',2);
ht=text(-.03,-.0,'1000 neurons','fontangle','normal','fontsize',8,'HorizontalAlignment','left');
ht.Rotation = 90;
axis off;
axis tight;

% covariance computation schematic
i=i+1;
hs{i}=my_subplot(5,3,10,[.83 .83]);
hs{i}.Position(1) = hs{i}.Position(1)-.03;
hs{i}.Position(2) = hs{i}.Position(2)-.05;
im = imread('cov_schematic2.png');
image(im);
axis off;
axis image;

% ----- PEER PRED PC TRACES ---------%
i=i+1;
hs{i}=my_subplot(5,3,11,[.75 .75]);
hs{i}.Position(1) = hs{i}.Position(1)-.05;
hs{i}.Position(2) = hs{i}.Position(2)-.06;
v2 = exampleV2 - mean(exampleV1,2);
v1 = exampleV1 - mean(exampleV1,2);
v2 = v2 ./ std(exampleV1(:,:),1,2);
v1 = v1 ./ std(exampleV1(:,:),1,2);
sk = skewness(v1,1,2);
v1 = v1 .* sign(sk);
v2 = v2 .* sign(sk);
t=0;
tr = 0+[1:75];
for n = [1 10 100 1000]
    t=t+1;
    plot(v1(n,tr)-(t-1)*6,'color',cpc(4*(t-1)+1,:));
    hold all;
    plot(v2(n,tr)-(t-1)*6,'color',.8*[1 1 1]);
    text(-0,.9-(t-1)*.22,sprintf('%d',n),'color',cpc(4*(t-1)+1,:),'fontangle','normal','HorizontalAlignment','right','fontsize',8);
end
text(.9,1.05,'test','color',0.8*[1 1 1],'fontangle','normal','HorizontalAlignment','right','fontsize',8);
text(.3,1.05,'train','color',cpc(1,:),'fontangle','normal','HorizontalAlignment','right','fontsize',8);
box off;
axis off;
text(0,1.05,'PC #','fontangle','normal','HorizontalAlignment','right','fontsize',8);
%set(gca,'xtick',2.^[0:4:10],'ytick',[0:.10:.30]);

% ----- PEER PRED PC ---------%
t=0;
for n = [1 10 100 1000]
	i=i+1;
	t=t+1;
	hs{i}=my_subplot(5,8,28+t,[xh*.9 yh*.9]);
	hs{i}.Position(2)=hs{i}.Position(2)-.07;
	if t==1
		hs{i}.Position(1)=hs{i}.Position(1)+.1;
	else
		hs{i}.Position(1)=hs{i}.Position(1)+.03*(4-t);
	end
	plot(exampleV1(n,:)*-1,exampleV2(n,:)*-1,'.','color',cpc(4*(t-1)+1,:),'markersize',2);
	box off;
	axis square;
	axis tight;
	if t==1
		xlabel('train');
		ylabel({'test'});
	end
	title(sprintf('PC %d \n r=%1.2f',n,corr(exampleV1(n,:)',exampleV2(n,:)')),'FontWeight','normal')
	if n==1
		set(gca,'xtick',[0 2000],'ytick',[0 2000]);
	else
		set(gca,'xtick',[0 2000],'ytick',[0 2000],'xticklabel',{},'yticklabel',{});
	end
	axis([-1000 2000 -1000 2000]);
end

% ----- power-law ---------%
i=i+1;
hs{i}=my_subplot(5,5,21,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.02;
hs{i}.Position(1)=hs{i}.Position(1)+.025;
for d = 1:ndat
	loglog(cov_neur(:,d)/sum(cov_neur(:,d)),'color',cdat(d,:),'linewidth',0.5)
	hold all;
end
loglog(nanmean(cov_neur,2)/sum(nanmean(cov_neur,2)),'k','linewidth',1)
box off;
axis square;
xlabel('PC dimension');
ylabel({'reliable','variance'});
set(gca,'xtick',10.^[0:3],'ytick',10.^[-4:2:10]);
axis tight;
ylim([1e-5 1]);
text(0.25,0.95,'\alpha = 1.14','fontsize',8)
grid on;
grid minor;
grid minor;

% ----- PEER PRED PC ---------%
i=i+1;
hs{i}=my_subplot(5,5,22,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.02;
hs{i}.Position(1)=hs{i}.Position(1)+.02;
for d = 1:ndat
    semilogx(cumsum(cov_neur(:,d))./sum(cov_neur(:,d)),'color',cdat(d,:),'linewidth',0.5);
    hold all;
end
semilogx(mean(cumsum(cov_neur)./sum(cov_neur),2),'k','linewidth',2);
box off;
axis square;
xlabel('PC dimension');
ylabel({'reliable','variance'});
title({'cumulative'},'fontweight','normal')
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]);
axis tight;
ylim([0 1]);

% ----- PEER PRED PC ---------%
i=i+1;
hs{i}=my_subplot(5,5,23,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.02;
hs{i}.Position(1)=hs{i}.Position(1)+.01;
for d = 1:ndat
    semilogx(cov_neur(:,d)./var_neur(:,d),'color',cdat(d,:),'linewidth',0.5);
    hold all;
end
semilogx(mean(cov_neur./var_neur,2),'color',.8*[1 1 1],'linewidth',2);
%expccum = cumsum(nanmean(cov_neur,2))./cumsum(nanmean(var_neur,2));
%semilogx(expccum,'k','LineWidth',2);
%text(0.2,1.15,'cumulative','FontSize',8);
text(0.05,0.25,'mean','FontSize',8,'color',.6*[1 1 1]);
box off;
axis square;
xlabel('PC dimension');
ylabel({'variance','explained'});
%title({'explainable','variance'},'fontweight','normal')
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]);
axis tight;
ylim([0 1]);



% ----- PC By PC VAR EXP ---------%
i=i+1;
hs{i}=my_subplot(5,5,24,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.02;
hs{i}.Position(1)=hs{i}.Position(1)+.01;
cm(4,:) = [0 0 0];
idim = [1 1 1 3];
ij = [1:3 7];
for j=1:4
    semilogx((nanmean(cov_neur,2) - nanmean(cov_res_beh(:,idim(j),:,ij(j)),3))./nanmean(var_neur,2),'color',cm(j,:),'linewidth',0.5);
    hold all;
end
semilogx(mean(cov_neur./var_neur,2),'color',.8*[1 1 1],'linewidth',2);
box off;
axis square;
xlabel('PC dimension');
ylabel({'variance','explained'});
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]);
axis tight;
axis([0 1024 0 1]);


% ----- PC By PC VAR EXP ---------%
i=i+1;
hs{i}=my_subplot(5,5,25,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.02;
hs{i}.Position(1)=hs{i}.Position(1)+.01;
cm(4,:) = [0 0 0];
idim = [1 1 1 3];
ij = [1:3 7];
exbeh = [];
for j=1:4
    exb = cumsum(cov_neur - squeeze(cov_res_beh(:,idim(j),:,ij(j))))./cumsum(var_neur);
	exbeh(j,:) = exb(128,:);
end
for d= 1:ndat
	plot(exbeh(:,d),'color',cdat(d,:));
	hold all;
end
plot(nanmean(exbeh,2),'k','linewidth',2);
%expccum = cumsum(nanmean(cov_neur,2))./cumsum(nanmean(var_neur,2));
%plot(ones(4,1)*expccum(128),'color',.8*[1 1 1])
box off;
axis square;
ylabel({'variance','explained'});
grid on;
grid minor;
grid minor;
set(gca,'xtick',[1:4],'ytick',[0:.05:1],'xticklabel',...
	{'\color[rgb]{.2,.8,.2} run','\color[rgb]{.5,.6,.5} pupil','\color[rgb]{.1,.9,.6}whisk','all 3'});
set(gca,'XTickLabelRotation',45)
axis tight;
title('PC 1-128','fontweight','normal');
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

%%
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


%%


% ----- PEER PRED TOTAL ---------%
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




% ------ PEER PRED DIAGRAM ------ %
i=i+1;
hs{i}=my_subplot(5,4,17,[xh yh]);
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

