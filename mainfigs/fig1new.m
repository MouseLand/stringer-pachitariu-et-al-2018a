function fig1new(matroot)

pc = load(fullfile(matroot,'corr1stpc.mat'));
clust = load(fullfile(matroot,'clust1D.mat'));

load(fullfile(matroot,'expv_behavior_PC.mat'));
load(fullfile(matroot,'PCpredict.mat'));

dex = 2;

%%
close all;
default_figure([1 1 4.75 7]);

%%
%trange = 750+[1:1500];
trange = [1:750];

pcplot=(my_conv2(zscore(pc.results.spks(1:2:end,trange),1,2),[3 0.5],[1 2]));
clustplot=(my_conv2(zscore(clust.results.spks(end:-1:1,trange),1,2),[3 0.5],[1 2]));
%%
clear hs;
i = 0;
cgray = colormap('gray');

dex = 2;

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
hs{i}.Position(1) = hs{i}.Position(1)+.02;
hs{i}.Position(2) = hs{i}.Position(2) + .01;
axis off;
hp=hs{i}.Position;
axes('position',[hp(1)-.08 hp(2)-.03 1.4*hp(3:4)]);
im = imread('planes.PNG');
imagesc(im);
title({'','','  imaging','setup'},'fontweight','normal')
set(gca,'fontsize',6);
axis off;
axis image;

i=i+1;
hs{i}=my_subplot(5,5,2,[xh yh]);
hs{i}.Position(1) = hs{i}.Position(1)+.0;
hs{i}.Position(2) = hs{i}.Position(2) + .01;
axis off;
hp=hs{i}.Position;
axes('position',[hp(1)-.07 hp(2)-.03 1.4*hp(3:4)]);
im = imread('explane.PNG');
imagesc(im);
title('plane 5/11','fontweight','normal','fontsize',6);
axis off;
axis image;

% ---------- correlation histogram -------------- %
i = i+1;
hs{i} = my_subplot(5,5,3,[xh yh]);
hs{i}.Position(1) = hs{i}.Position(1)-.015;
hs{i}.Position(2) = hs{i}.Position(2) + .01;
bins = [-.3:.001:.75];
nb = histcounts(pc.results.cshuff(1:1:end),bins);
%nb = nb/sum(nb);
%area(bins(1:end-1),log(nb),-17,'facecolor',.5*[1 1 1],'edgecolor','none','showbaseline','off')
pbins = bins(1:end-1)+(bins(2)-bins(1))/2;
plot(pbins, log10(nb), 'color', .5*[1 1 1]);
hold all;
nb = histcounts(pc.results.ccall(1:1:end),bins);
%nb = nb/sum(nb);
plot(pbins,log10(nb),'color',cdat(dex,:))
text(.4,.9,'shuffled','fontangle','normal','fontsize',6,'color',.5*[1 1 1]);%,'fontweight','bold');
text(.4,1.1,'original','fontangle','normal','fontsize',6,'color',cdat(dex,:));%,'fontweight','bold');
axis tight;
%ylim([-16 -3]);
xlim([-.25 .7]);
set(gca,'ytick',[0:2:6],'yticklabel',{'10^0','10^2','10^4',...
    '10^6'});
set(gca,'fontsize',6);
xlabel('correlation');
ylabel('# of pairs');
box off;
axis square;

% ---------- mean and std of correlations -------------- %
i = i+1;
hs{i} = my_subplot(5,5,4,[xh yh]);
hs{i}.Position(1) = hs{i}.Position(1)+.0;
hs{i}.Position(2) = hs{i}.Position(2) + .01;
hold all;
id = [2 1 5 3 4 6:9]; % show gray screen first
for d = 1:ndat
    errorbar(d, pc.results.mcorr(id(d)), pc.results.stdcorr(id(d),1), pc.results.stdcorr(id(d),2),'.', 'color', cdat(id(d),:))
end
hi = .23;
y1 = [.5 3.7];
y2 = [3.3 9.5];
%for j = 1:2
%plot([y1(j) y2(j)],hi*[1 1],'k');
%plot(y1(j)*[1 1],[hi-.03 hi],'k');
%plot(y2(j)*[1 1],[hi-.03 hi],'k');
%end
text(.1,1.4,{' gray','screen'},'fontsize',6,'fontangle','normal','color',cdat(2,:))
text(0.7,1.4,{'dark-','ness'},'fontsize',6,'fontangle','normal','color',cdat(9,:))
ylabel('correlation');
xlabel('recordings');
set(gca,'fontsize',6);
set(gca,'xtick',[]);
box off;
axis square;
ylim([-.05 .25]);

pcs = pc.results.firstpcs(:,:);
pcs = zscore(pcs.*sign(skewness(pcs)),1,1);
behall = zscore(pc.results.behavior(:,:));
behall(:,3) = behall(:,3)*-1;
behall = behall - min(behall,[],1);

% ----------- correlation with behavior ----------%
tstr = {'running','pupil area','whisking','all 3'};
behall0 = min(5,behall);

i = i+1;
hs{i} = my_subplot(5,5,5,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + .01;
hs{i}.Position(1) = hs{i}.Position(1)+.015;
%hs{i}.Position(2) = hs{i}.Position(2)+0.05;
hold all; 
plot(pcs(:,1),behall(:,1),'.','color',cdat(dex,:),'markersize',4);
xlabel('magnitude');
ylabel({'running','speed (a.u.)'});
axis square;
axis tight;
title({'PC 1','r=0.75'},'fontweight','normal','fontsize',6);
set(gca,'fontsize',6);
%axis([0 6 0 6]);
ylim([0 7]);


% --------- PC img ------------------ %
mpcs = min(pcs,[],1);
%pcs = pcs - mpcs;
i = i+2;
x0 = .26;
y0 = .76;
xh0 = .7;
yh0 = .4;
hs{i} = my_subplot(5,1,2,[.92 .7]);%axes('position',[x0 y0 xh0 yh0]);
hs{i}.Position(2)=hs{i}.Position(2)-0.09;
hs{i}.Position(1)=hs{i}.Position(1);
%hs{i}.Position(1)=hs{i-4}.Position(1);
pos=hs{i}.Position;
[NN,NT] = size(pcplot);
imagesc(pcplot(:,:), [0 .3]);
hold all;
set(gca,'fontsize',6);
colormap(hs{i},flipud(cgray));
plot(-5*[1 1], [NN-500 NN],'k','linewidth',2);
ht=text(-.03,-.0,'1000 neurons','fontangle','normal','fontsize',6,'HorizontalAlignment','left');
ht.Rotation = 90;
%ht=text(-0.03,.5,'sorted by 1st PC weights','fontangle','normal','fontsize',6,...
%    'HorizontalAlignment','center','fontweight','bold');
%ht.Rotation = 90;
axis off;
axis tight;

tb=1;
i=i-1;
hs{i}=axes('position',[pos(1) pos(2)+pos(4)+.03 pos(3)*tb .11]);
hold all;
cm(1,:) = [.2 .8 .2];
cm(2,:) = [.5 .6 .5];
cm(3,:) = [.1 .9 .6];
tstr = {'running','pupil area','whisking'};
for j = 1:3
	beh = behall0(trange(1:length(trange)/3*tb),j);
	plot(beh,'linewidth',1,'color',cm(j,:),'linewidth',1);
	text(.03,1.-(j-1)*0.12,tstr{j},'fontangle','normal','fontsize',6,'color',cm(j,:));
		%'fontweight','bold');
end
plot(pcs(trange(1:length(trange)/3*tb),1),'-.','color',cpc(1,:),'linewidth',0.5);
text(.4,.95,'PC 1','HorizontalAlignment','right','fontangle','normal','color',cpc(1,:),'fontsize',6)
text(0,-.05,'10 s','fontangle','normal','fontsize',6);
plot([0 10/1.2],-1.5*[1 1],'k','LineWidth',2);
axis tight;
set(gca,'XColor','w');
%xlim([-5/3 length(trange)/3]);


% ---------- CLUST img ------------- %
i=i+2;
hs{i} = my_subplot(5,1,3,[1 0.9]);%axes('position',[x0 y0 xh0 yh0]);
hs{i}.Position(2)=hs{i}.Position(2)-.03;
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
ht=text(-.03,-.0,'1000 neurons','fontangle','normal','fontsize',6,'HorizontalAlignment','left');
ht.Rotation = 90;
axis off;
axis tight;

% covariance computation schematic
i=i+1;
hs{i}=my_subplot(5,3,10,[.88 .88]);
hs{i}.Position(1) = hs{i}.Position(1)-.03;
hs{i}.Position(2) = hs{i}.Position(2)-.02;
im = imread('cov_fig.png');
image(im);
axis off;
axis image;

% ----- PEER PRED PC TRACES ---------%
i=i+1;
hs{i}=my_subplot(5,3,11,[.82 .75]);
hs{i}.Position(1) = hs{i}.Position(1)-.01;
hs{i}.Position(2) = hs{i}.Position(2)-.04;
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
    text(-0,.9-(t-1)*.22,sprintf('%d',n),'color',cpc(4*(t-1)+1,:),'fontangle','normal','HorizontalAlignment','right','fontsize',6);
end
text(1,1.1,'cell set 2','color',0.8*[1 1 1],'fontangle','normal','HorizontalAlignment','right','fontsize',6);
text(.4,1.1,'cell set 1','color',cpc(1,:),'fontangle','normal','HorizontalAlignment','right','fontsize',6);
plot([0 10/1.2],[1 1]-22,'k')
text(0,0,'10 s','fontangle','normal','fontsize',6);
box off;
axis off;
axis tight;
text(0,1.0,'SVC #','fontangle','normal','HorizontalAlignment','right','fontsize',6);
set(gca,'fontsize',6);
%set(gca,'xtick',2.^[0:4:10],'ytick',[0:.10:.30]);

% ----- PEER PRED PC ---------%
t=0;
i=i+1;
clear htpos;
	
for n = [1 10 100 1000]
	t=t+1;
	if t==1
		hs{i}=my_subplot(9,5,34+mod(t-1,2)+(t>2)*5,[xh*.9 yh*.9]);
		hs{i}.Position(2)=hs{i}.Position(2)+.02;
		hs{i}.Position(1)=hs{i}.Position(1)+.03;
	else
		hp=my_subplot(9,5,34+mod(t-1,2)+(t>2)*5,[xh*.9 yh*.9]);
		hp.Position(2)=hp.Position(2)+.02;
		hp.Position(1)=hp.Position(1)-.04*mod(t-1,2)+.03;
	end
	
	plot(exampleV1(n,:)*-1,exampleV2(n,:)*-1,'.','color',cpc(4*(t-1)+1,:),'markersize',2);
	box off;
	axis square;
	axis tight;
	if t==1
		xlabel('  cell set 1');
		ylabel({'  cell set 2'});
	end
	text(0.1,1.7,sprintf('SVC %d \n r=%1.2f',n,corr(exampleV1(n,:)',exampleV2(n,:)')),'FontAngle','normal','fontsize',6);
	%if n==1
	%	set(gca,'xtick',[0 2000],'ytick',[0 2000]);
	if 1
		set(gca,'xtick',[0 2000],'ytick',[0 2000],'xticklabel',{},'yticklabel',{});
	end
	axis([-1000 2000 -1000 2000]);
	set(gca,'fontsize',6);
end

% ----- PEER PRED PC ---------%
i=i+1;
hs{i}=my_subplot(5,4,22-5,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.03;
hs{i}.Position(1)=hs{i}.Position(1)+.02;
for d = 1:ndat
    semilogx(cov_neur(:,d)./var_neur(:,d) * 100,'color',cdat(d,:),'linewidth',0.5);
    hold all;
end
semilogx(mean(cov_neur./var_neur,2) * 100,'color',.8*[1 1 1],'linewidth',2);
%expccum = cumsum(nanmean(cov_neur,2))./cumsum(nanmean(var_neur,2));
%semilogx(expccum,'k','LineWidth',2);
%text(0.2,1.15,'cumulative','FontSize',8);
text(0.5,1,'mean','FontSize',6,'color',.6*[1 1 1]);
box off;
axis square;
xlabel('SVC dimension');
ylabel('% reliable variance');
%title({'explainable','variance'},'fontweight','normal')
set(gca,'xtick',10.^[0:2:3],'ytick',[0:50:100]);
set(gca,'fontsize',6);
axis tight;
ylim([0 100]);

% ----- power-law ---------%
i=i+1;
hs{i}=my_subplot(5,4,23-5,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.03;
hs{i}.Position(1)=hs{i}.Position(1)+.02;
for d = 1:ndat
	loglog(cov_neur(:,d)/sum(cov_neur(:,d)),'color',cdat(d,:),'linewidth',0.5)
	hold all;
end
loglog(nanmean(cov_neur,2)/sum(nanmean(cov_neur,2)),'k','linewidth',1)
box off;
axis square;
xlabel('SVC dimension');
ylabel('reliable variance');
set(gca,'xtick',10.^[0:2:3],'ytick',10.^[-4:2:10]);
set(gca,'fontsize',6);
axis tight;
ylim([1e-5 1]);
text(0.15,1.03,'\alpha = 1.14','fontsize',6)
grid on;
grid minor;
grid minor;

% ----- PEER PRED PC ---------%
% i=i+1;
% hs{i}=my_subplot(5,5,23,[xh yh]);
% hs{i}.Position(2)=hs{i}.Position(2)-.03;
% hs{i}.Position(1)=hs{i}.Position(1)+.01;
% for d = 1:ndat
%     semilogx(cumsum(cov_neur(:,d))./sum(cov_neur(:,d)),'color',cdat(d,:),'linewidth',0.5);
%     hold all;
% end
% semilogx(mean(cumsum(cov_neur)./sum(cov_neur),2),'k','linewidth',2);
% box off;
% axis square;
% xlabel('PC dimension');
% ylabel({'reliable variance'});
% title({'cumulative'},'fontweight','normal')
% set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]);
% axis tight;
% ylim([0 1]);

% ----- PC By PC VAR EXP ---------%
i=i+1;
hs{i}=my_subplot(5,4,24-5,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.03;
hs{i}.Position(1)=hs{i}.Position(1)+.02;
cm(4,:) = [0 0 0];
idim = [1 1 1 3];
ij = [1:3 7];
for j=1:4
    semilogx(100*(nanmean(cov_neur,2) - nanmean(cov_res_beh(:,idim(j),:,ij(j)),3))./nanmean(var_neur,2),'color',cm(j,:),'linewidth',0.5);
    hold all;
end
semilogx(100*mean(cov_neur./var_neur,2),'color',.8*[1 1 1],'linewidth',2);
box off;
axis square;
xlabel('SVC dimension');
ylabel('% reliable variance');
text(.08,.55,'run','fontangle','normal','fontsize',6,'color',cm(1,:));
text(.13,.4,'pupil','fontangle','normal','fontsize',6,'color',cm(2,:));
text(.23,.25,'whisk','fontangle','normal','fontsize',6,'color',cm(3,:));
text(.03,.7,'all 3','fontangle','normal','fontsize',6,'color',cm(4,:));
text(0.5,1,'mean','FontSize',6,'color',.6*[1 1 1]);
%ylabel({'variance','explained'});
%ht=text(-.85,-.45,{'% reliable var'},'fontsize',6,'fontangle','normal');
%ht.set('rotation',90);
set(gca,'xtick',10.^[0:2:3],'ytick',[0:50:100]);
set(gca,'fontsize',6);
axis tight;
axis([0 1024 0 100]);


% ----- PC By PC VAR EXP ---------%
i=i+1;
hs{i}=my_subplot(5,4,25-5,[xh yh]);
hs{i}.Position(2)=hs{i}.Position(2)-.03;
hs{i}.Position(1)=hs{i}.Position(1)+.02;
cm(4,:) = [0 0 0];
idim = [1 1 1 3];
ij = [1:3 7];
exbeh = [];
for j=1:4
    exb = cumsum(cov_neur - squeeze(cov_res_beh(:,idim(j),:,ij(j))))./cumsum(var_neur);
	exbeh(j,:) = exb(128,:);
end
for d= 1:ndat
	plot(100* exbeh(:,d),'color',cdat(d,:));
	hold all;
end
plot(100*nanmean(exbeh,2),'k','linewidth',2);
%expccum = cumsum(nanmean(cov_neur,2))./cumsum(nanmean(var_neur,2));
%plot(ones(4,1)*expccum(128),'color',.8*[1 1 1])
box off;
axis square;
ylabel('% reliable variance');
%ylabel({'variance','explained'});
grid on;
grid minor;
grid minor;
set(gca,'xtick',[1:4],'ytick',[0:5:15],'xticklabel',...
	{'\color[rgb]{.2,.8,.2} run','\color[rgb]{.5,.6,.5} pupil','\color[rgb]{.1,.9,.6}whisk','all 3'});
set(gca,'XTickLabelRotation',45)
axis tight;
ylim([0 15])
title('SVC 1-128','fontweight','normal','fontname','Helvetica');
set(gca,'fontsize',6);
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
%     text(.03,1-.13*(j-1),tstr{j},'color',cm(j,:),'fontangle','normal','fontsize',6);
%     xlim([0 1]);
%     axis square;
%     axis tight;
% end

tstr{1} = '';
for t = 2:30
	tstr{t} = '';
end

% -------------- LETTERS
hp=.08;
hy=1.27;
deffont=8;
for j = [1:length(hs)]
	hp0=hp; hy0=hy;
    if j==6
        hp0 = 0.04;
        hy0 = 1.15;
	elseif j==8 || j==7
		hy0 = 1.02;
		hp0 = .04;
	elseif j==9
		hp0=-.01;
		hy0=0.95;
	elseif j==10
		hp0=0.07;
		hy0=1.16;
	elseif j==11
		hp0=0.05;
		hy0=1.7;		
    elseif j>=12
		hy0=hy*.98;
		hp0=.1;
	else
        hp0=hp;
        hy0 = hy;
		
    end
    hpos = hs{j}.Position;
    lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
    axes('position', lpos);
    t=text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','FontName', 'Helvetica');
    axis([0 1 0 1]);
    axis off;
        
    if j==10
        tsh=.02;
    else
        tsh=.03;
    end
    axes('position', [lpos(1)+tsh lpos(2)-.006 lpos(3) lpos(4)]);
    text(0,0, tstr{j},'fontangle','normal','fontsize',6);
    axis([0 1 0 1]);
    axis off;
    
end

%%

print(fullfile(matroot,'fig1new.pdf'),'-dpdf');



