function fig2newnew(matroot)

load('exampleMovie.mat');
load(fullfile(matroot,'example_behavior.mat'));
load(fullfile(matroot,'PCpredict.mat'));
load(fullfile(matroot,'expv_behavior_PC.mat'));
try
    load(fullfile(matroot,'expv_timedelay.mat'));
catch
end

%%
close all;
figX = 4.75;
figY = 4;
default_figure([1 1 figX figY]);
%%

yh = .05;
xh = .02;

dex=2;
clf;
clear hs;
cm = colormap('winter');
load('cdat.mat');
% cd=colormap('hot');
% clear cdat;
% cdat(1,:) = [.9 .8 .6];
% cdat(2,:) = [.9 .9 0];
% cdat(5,:) = [0.8 .7 .3];
% cred=cd(12+ [1:4:4*6],:);
% rng(1);
% cdat([3 4 6:9],:) = cred(randperm(6),:);
%cdat(3,:) = cdat(4,:);
%load('cdat.mat');
ifr =1;
cgray=colormap('gray');

ncomps = 4;
cm = cm(round(linspace(1,size(cm,1),ncomps)),:);

nPCs = 4;
cpc = max(0,colormap('spring'));
cpc = cpc(1:32,:);
cpc = cpc(round(linspace(1,size(cpc,1),nPCs)),:);
spc = [-1 1 -1 -1];
ndat = size(cdat,1);

i = 0;

i=i+1;
[ny,nx,~]=size(exframes);
hs{i} = my_subplot(4,5,1,[.65 .6]);
hs{i}.Position(1) = hs{i}.Position(1)+.01;
hs{i}.Position(2) = hs{i}.Position(2)+.01;
hold all;
imagesc([1:nx]-50,[1:ny]-50,exframes(:,:,ifr));
imagesc([1:nx],[1:ny],exframes(:,:,ifr+1));
text(-.1,0.22,'t','fontsize',6);
text(0,0,'t+1','fontsize',6);
colormap(gca,cgray);
set(gca,'ydir','reverse');
axis image;
axis off;

i=i+1;
hs{i} = my_subplot(4,5,2,[.65 .6]);
hs{i}.Position(1) = hs{i}.Position(1)+.03;
imagesc(abs(exframes(:,:,ifr+1)-exframes(:,:,ifr)),[-150 150]);
colormap(gca,redblue);
axis image;
axis off;
text(0.03,1.24,'motion energy','fontsize',6);

x0 = hs{1}.Position(1) + hs{1}.Position(3);
x1 = hs{2}.Position(1);
y0 = hs{2}.Position(2) + hs{2}.Position(4)/2;
annotation('arrow',[x0 x1],[y0 y0]);


[ny,nx,~]=size(exmotionMask);
for j = 1:ncomps-1
    hp = my_subplot(4,5,2+j,[.65 .6]);
    hp.Position(1) = hp.Position(1)-.055*(j-1) + .05;
    if j==1
        i=i+1;
        hs{i}=hp;
    end
    %axes('position',[xc yc+yh*(ncomps-j) .1 .12]);
    pmask = exmotionMask(:,:,j);
    pmask = pmask * spc(j);
    hold all;
    text(.4,1.24,sprintf('PC %d',j),'color',cm(j,:),'fontsize',6);
    imagesc(pmask,20e-3*[-1 1]);
    patch([0 0 nx+1 nx+1], [0 ny+1 ny+1 0],'k','edgecolor',cm(j,:),...
        'linewidth',1,'facecolor','none');
    colormap(hp,redblue);
    set(gca,'ydir','reverse');
    axis off;
    axis image;
end
axes('position',[hp.Position(1)+hp.Position(3)-.02 hp.Position(2)+.02 .02 hp.Position(4)-.04]);
cb=colorbar('location','eastoutside');
cb.Ticks=[0 .5 1];
cb.TickLabels={'-0.2','0','0.2'};
axis off;
axis tight;
colormap(gca,'redblue');
    

x0 = hs{2}.Position(1) + hs{2}.Position(3);
x1 = hs{3}.Position(1);
y0 = hs{3}.Position(2) + hs{3}.Position(4)/2;
annotation('arrow',[x0 x1],[y0 y0]);

% ---------- SCHEMATIC ------------------%
i=i+1;
hs{i} = my_subplot(4,5,6,[1 .6]);
hs{i}.Position(2)=hs{i}.Position(2)+.04;
hs{i}.Position(1) = hs{i}.Position(1) + .04;
hold all;
trange = 200+[1:100];
NT = length(trange);
for j = 1:ncomps
    xp=zscore(xtest{dex}(j,trange));
    xp = xp * spc(j);
    plot(xp-j*8, 'color', cm(j,:));
end
plot([0 30/1.2],[1 1]*(-j*8-2),'k');
text(0,0.0,' 30 s','fontsize',6,'fontangle','normal');
axis tight;
axis off;
text(0,1.2,'motion energy PCs','fontsize',6);
text(1.95,1.2,'neural SVCs','fontsize',6);
text(1.35,1,'rank','fontsize',6);
h2=my_subplot(4,5,8,[1 .6]);
h2.Position(1) = h2.Position(1) +.04;
h2.Position(2) = hs{i}.Position(2) + .0;
hold all;
cp = .4*[1 1 1];
ipc = [1 3 4 5];
for j = 1:nPCs
    xp=vtest{dex}(ipc(j),trange);
    pred=vpred{dex}(ipc(j),trange);
    pred = pred / std(xp);
    xp = xp / std(xp);
	ssk = sign(skewness(xp));
	xp = xp * ssk;
	pred = pred*ssk;
    r = corr(vtest{dex}(j,:)',vpred{dex}(j,:)');
    plot(NT*3+[1:NT],+xp-j*8, 'color', cpc(j,:));
    plot(NT*3+[1:NT],pred-j*8, 'color', cp);
    %text(NT*3.6,-j*8+8,sprintf('r^2=%1.2f',r^2),'units','data','fontsize',6,...
    %    'fontweight','bold','color',cpc(j,:));
end
axis tight;
axis off;

% ---- BALLS ----%
nHidden = 2;
call{1} = cm;
call{2} = repmat(cp,nHidden,1);
call{3} = cpc;
pos(1) = hs{i}.Position(1)+hs{i}.Position(3)*1.1;
pos(3) = h2.Position(1)*.95 - pos(1);
pos(2) = hs{i}.Position(2)+.01;
pos(4) = hs{i}.Position(4)*.88;
nls = [ncomps nHidden nPCs];
msize = 6;
layeredNet(pos, nls, call, msize,[0 0 0]);
% little ellipses
axes('position',[pos(1) pos(2)-.05 pos(3) .06]);
ipos = [1 2 3];
hold all;
for j = 1:3
    if j ==2
        jrange=[4:6];
    else
        jrange=[1:3];
    end
    plot(ones(1,3)*ipos(j),jrange,'k.');
end
axis tight;
axis off;


xh=.65;
yh=.65;

% ------------ PC VAR EXP ------------%
i=i+1;
hs{i} = my_subplot(4,4,8,[xh*1.1 yh*1.1]);
hs{i}.Position(1) = hs{i}.Position(1) -.06;
hs{i}.Position(2) = hs{i}.Position(2)+.06;
hold all;
shadedErrorBar([1:1024],100*mean((cov_neur-squeeze(cov_res_beh(:,6,:,9)))./var_neur,2),...
	std(100*(cov_neur-squeeze(cov_res_beh(:,6,:,9)))./var_neur,1,2)/sqrt(ndat-1),'lineProps',{'color','b','linewidth',1});
shadedErrorBar([1:1024],100*mean(cov_neur./var_neur,2),...
	std(100*cov_neur./var_neur,1,2)/sqrt(ndat-1),'lineProps',{'color',.8*[1 1 1],'linewidth',1});
shadedErrorBar([1:1024],100*mean((cov_neur-squeeze(cov_res_beh(:,3,:,7)))./var_neur,2),...
	std(100*(cov_neur-squeeze(cov_res_beh(:,3,:,7)))./var_neur,1,2)/sqrt(ndat-1),'lineProps',{'color',[0 0.2 0],'linewidth',1});
%axis([0 .06 0 .06]);
set(gca,'xtick',10.^[0 1 2])
text(.55,.9,{'max explainable'},'color',.6*[1 1 1],'fontsize',6)
text(0.5,0.25,{'run+pupil+whisk'},'fontangle','normal','fontsize',6,'color',[0 .2 0]);
text(.3,0.6,'face','fontangle','normal','fontsize',6,'color','b');
text(0.5,-0.25,'SVC dimension','HorizontalAlignment','center');
ht=text(-.4,0.5,{'% variance','explained'},'HorizontalAlignment','center','VerticalAlignment','middle');
set(ht,'rotation',90);
set(gca,'fontsize',6);
box off;
axis square;
axis tight;
ylim([0 100]);
xlim([0 128]);
set(gca,'xscale','log');
grid on;
grid minor;
grid minor;


% ------------ RASTERS
trange = 190+[1:350];%[800:1300];
%trange = [trange(5:end-13) 830+[1:45]];% 1200+[1:50]];
i=i+1;
hs{i}=my_subplot(4,2,5,[.95 .8]);
hs{i}.Position(2) = hs{i}.Position(2) +.01; 
hs{i}.Position(1) = hs{i}.Position(1) +.02; 
yt = ytest{dex}(isortembed{dex}(end:-1:1),trange);
ytstd = max(1e-3,std(yt,1,2));
yt = yt ./ ytstd;
imagesc(my_conv2(yt, 6, 1),[0 .3])
colormap(hs{i},flipud(cgray));
axis off;
text(0,1.15,'Neural activity (test set) sorted by 1D embedding','fontangle','normal','fontsize',6);
NN=size(yt,1);
hold all;
plot([-2 10/1.2], (NN+100)*[1 1],'k');
plot([0 0]-2, NN+[0 -1000],'k');
ht=text(-.05,0,'1000 neurons','fontsize',6,'fontangle','normal');
ht.Rotation=90;
text(0,0,'10 s','fontsize',6,'fontangle','normal');
axis tight;


hp=my_subplot(4,2,7,[.95 .8]);
%hp.Position(1) = hs{i}.Position(1);
%hp.Position(3) = hs{i}.Position(3);
%hp.Position(1) = hs{4}.Position(1);
hp.Position(1) = hs{i}.Position(1);
%hp.Position(4) = .22;
yp = ypred{dex}(isortembed{dex}(end:-1:1),trange);
yp = yp ./ ytstd; % scale prediction by same amount
yp = my_conv2(yp, 0.5, 1);
imagesc(yp,[0 .3])
colormap(hp,flipud(cgray));
axis off;
text(0,1.15,'Neural activity prediction (test set) from face motion','fontangle','normal','fontsize',6);
ht=text(-.05,0,'1000 neurons','fontsize',6,'fontangle','normal');
ht.Rotation=90;
hold all;
plot([-2 10/1.2], (NN+100)*[1 1],'k');
plot([0 0]-2, NN+[0 -1000],'k');
text(0,0,'10 s','fontsize',6,'fontangle','normal');
axis tight;



% --------------------- DIMS VS VAR EXP -------------- %
i=i+1;
hs{i} = my_subplot(4,5,14,[xh yh]);
hs{i}.Position(1) = hs{i}.Position(1) - .035;
hs{i}.Position(2) = hs{i}.Position(2)+.03;
%hs{i}.Position(2) = hs{i-1}.Position(2);
expcbeh = cumsum(cov_neur - permute(squeeze(cov_res_beh(:,:,:,9)),[1 3 2]))./cumsum(var_neur);
expcbeh = squeeze(expcbeh(128,:,:))';
for j = 1:ndat
    semilogx(ndims0,100*expcbeh(:,j),'color',cdat(j,:));
    hold all;
end
box off;
set(gca,'xtick',[1 4 16 64]);
text(0.5,-0.3,'rank','HorizontalAlignment','center');
ht=text(-.5,0.5,{'% variance','explained'},'HorizontalAlignment','center','VerticalAlignment','middle');
set(ht,'rotation',90);
axis tight;
ylim([0 30]);
xlim([0 128]);
hold all;
semilogx(ndims0,mean(expcbeh,2)*100,'b','linewidth',2);
expcbeh = cumsum(cov_neur - squeeze(cov_res_beh(:,3,:,7)))./cumsum(var_neur);
semilogx(ndims0,100*mean(expcbeh(128,:))*ones(length(ndims0),1),'color',[0 .2 0],'linewidth',1);
text(.1,.3,'run+pupil+whisk','fontangle','normal','fontsize',6,'color',[0 .2 0]);
axis square;
grid on;
grid minor;
grid minor;
title('SVC 1-128','fontweight','normal');
set(gca,'fontsize',6);

%-------------- PEER PRED -----------------%
i=i+1;
hs{i}=my_subplot(4,5,15,[xh yh]);
hs{i}.Position(2) = hs{i-1}.Position(2);
hs{i}.Position(1) = hs{i}.Position(1) - .0;
expcbeh = cumsum(cov_neur - squeeze(cov_res_beh(:,6,:,9)))./cumsum(var_neur);
expccum = cumsum(cov_neur)./cumsum(var_neur);
ev = 100* [expccum(128,:)' expcbeh(128,:)'];
hold all;
for j = 1:size(ev,1)
    plot(ev(j,1), ev(j,2), 'wo',...
        'markerfacecolor',cdat(j,:),'markersize',4);
end
plot([0 80],[0 80],'k','linewidth',1)
axis square;
box off;
ylim([0 80])
xlim([0 80]);
set(gca,'xtick',[0:20:80],'ytick',[0:20:80]);
text(0.5,-0.25,'max explainable','HorizontalAlignment','center');
ht=text(-.4,0.5,{'face'},'HorizontalAlignment','center','VerticalAlignment','middle');
set(ht,'rotation',90);
%ht=title({'% variance explained',''},'fontsize',6);
%ht.Position(2) = ht.Position(2)-4;
grid on;
grid minor;
grid minor;
%title('SVC 1-128','fontweight','normal');
set(gca,'fontsize',6);

% ------------ FACE (16D) + arousal ------------%
i=i+1;
hs{i} = my_subplot(4,5,19,[xh yh]);
hs{i}.Position(1) = hs{i-2}.Position(1);
hs{i}.Position(2) = hs{i}.Position(2)+.02;
hold all;
plot([0 30],[0 30],'k');
expcbeh3d = cumsum(cov_neur - squeeze(cov_res_beh(:,6,:,10)))./cumsum(var_neur);
for d = 1:ndat
    plot(100*expcbeh(128,d), 100*expcbeh3d(128,d),'wo',...
        'markerfacecolor',cdat(d,:),'markersize',4);
end
axis([0 .3 0 .3]*100);
text(0.5,-0.25,'face + 3D arousal','HorizontalAlignment','center');
ht=text(-.4,0.5,{'face'},'HorizontalAlignment','center','VerticalAlignment','middle');
set(ht,'rotation',90);box off;
axis square;
grid on;
grid minor;
grid minor;
%title('SVC 1-128','fontweight','normal');
set(gca,'fontsize',6);


% ---------- TIMELAG --------%
i=i+1;
hs{i} = my_subplot(4,5,20,[xh yh]);
hs{i}.Position(2) = hs{i-1}.Position(2);
evtlag = cumsum(tlag_cov_neur - tlag_cov_res_beh)./...
	cumsum(tlag_var_beh);
hold all;
for d = 1:ndat
	plot(tdelay,100*evtlag(128,:,d),'color',cdat(d,:),'linewidth',0.5)
end
plot(tdelay,mean(100*evtlag(128,:,:),3),'b','linewidth',2);
%axis tight;
text(0.5,-0.25,'time from behav. (s)','HorizontalAlignment','center');
ht=text(-.5,0.5,{'% variance','explained'},'HorizontalAlignment','center','VerticalAlignment','middle');
set(ht,'rotation',90);
%title('SVC 1-128','fontweight','normal');
set(gca,'fontsize',6);
grid on;
grid minor;
grid minor;



% -------------- LETTERS
hp=.065;
hy=1.15;
deffont=8;
for j = [1:length(hs)]
	if j==2 || j==3
		hp0=0.04;
		hy0=1.27;
	elseif j== 5
    	hp0=0.065;
		hy0=1.1;
    elseif j==4
        hp0 =0.03;
        hy0 = 1.3;
	elseif j==1
		hp0=0.04;
		hy0=1.2;
	%elseif j==5
    %    hy0=1.1;
    %    hp0=.065;
    elseif j==6
        hp0=.03;
        hy0=1.2;
	elseif j==7 || j==10
		hy0=hy;
		hp0=0.09;
		
	else
        hp0=hp;
        hy0 = hy;
    end
    hpos = hs{j}.Position;
    axes('position', [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01]);
    text(0,0, char(64+j),'fontsize',10,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
end

tl=squeeze(evtlag(128,:,:));

tli = interp1(tdelay,tl,[-8:.01:8]);

[~,ix1]=min(abs(max(tli,[],1)/2 - tli(1:800,:)));
[~,ix2]=min(abs(max(tli,[],1)/2 - tli(801:end,:)));

fwhm = (800-ix1+ix2)*.01

%%
print(fullfile(matroot,'fig2new.pdf'),'-dpdf');