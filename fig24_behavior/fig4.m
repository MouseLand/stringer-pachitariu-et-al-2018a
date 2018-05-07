function fig4(matroot)

load('exampleMovie.mat');
load(fullfile(matroot,'example_behavior.mat'));
load(fullfile(matroot,'expv_behavior_neurons.mat'));
try
    load(fullfile(matroot,'PCApred.mat'));
catch
    load('PCApred.mat');
end
try
    load(fullfile(matroot,'expv_timedelay.mat'));
catch
end

%%

figX = 6;
figY = 8;
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
ndat = size(expv_behavior,2);
cdat = cdat(round(linspace(1,size(cdat,1),size(expv_behavior,2))),:);

i = 0;

i=i+1;
[ny,nx,~]=size(exframes);
hs{i} = my_subplot(5,5,1,[.65 .6]);
hs{i}.Position(1) = hs{i}.Position(1)+.01;
hold all;
imagesc([1:nx]-50,[1:ny]-50,exframes(:,:,ifr));
imagesc([1:nx],[1:ny],exframes(:,:,ifr+1));
text(-.1,0.22,'t','fontsize',8);
text(0,0,'t+1','fontsize',8);
colormap(gca,cgray);
set(gca,'ydir','reverse');
axis image;
axis off;

i=i+1;
hs{i} = my_subplot(5,5,2,[.65 .6]);
hs{i}.Position(1) = hs{i}.Position(1)+.03;
imagesc(abs(exframes(:,:,ifr+1)-exframes(:,:,ifr)),[-150 150]);
colormap(gca,redblue);
axis image;
axis off;
text(0.03,1.24,'motion energy','fontsize',8);

x0 = hs{1}.Position(1) + hs{1}.Position(3);
x1 = hs{2}.Position(1);
y0 = hs{2}.Position(2) + hs{2}.Position(4)/2;
annotation('arrow',[x0 x1],[y0 y0]);


[ny,nx,~]=size(exmotionMask);
for j = 1:ncomps-1
    hp = my_subplot(5,5,2+j,[.65 .6]);
    hp.Position(1) = hp.Position(1)-.045*(j-1) + .05;
    if j==1
        i=i+1;
        hs{i}=hp;
    end
    %axes('position',[xc yc+yh*(ncomps-j) .1 .12]);
    pmask = exmotionMask(:,:,j);
    pmask = pmask * spc(j);
    hold all;
    text(.4,1.24,sprintf('PC %d',j),'color',cm(j,:),'fontsize',8);
    imagesc(pmask,20e-3*[-1 1]);
    patch([0 0 nx+1 nx+1], [0 ny+1 ny+1 0],'k','edgecolor',cm(j,:),...
        'linewidth',2,'facecolor','none');
    colormap(hp,redblue);
    set(gca,'ydir','reverse');
    axis off;
    axis image;
end
axes('position',[hp.Position(1)+hp.Position(3)-.005 hp.Position(2)+.02 .02 hp.Position(4)-.04]);
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
hs{i} = my_subplot(5,3,4,[.8 .6]);
hs{i}.Position(2)=hs{i}.Position(2)+.02;
hold all;
trange = 200+[1:100];
NT = length(trange);
for j = 1:ncomps
    xp=zscore(xtest{dex}(j,trange));
    xp = xp * spc(j);
    plot(xp-j*8, 'color', cm(j,:));
end
plot([0 30/1.2],[1 1]*(-j*8-2),'k');
text(0,0.0,' 30 s','fontsize',8,'fontangle','normal');
axis tight;
axis off;
text(0,1.2,'motion energy PCs','fontsize',8);
text(1.95,1.2,'neural PCs','fontsize',8);
text(3.05,1.2,'neurons','fontsize',8);

h2=my_subplot(5,3,6,[.8 .6]);
h2.Position(1) = h2.Position(1) - .14;
h2.Position(2) = h2.Position(2) + .02;
hold all;
cp = .4*[1 1 1];
for j = 1:nPCs
    xp=vtest{dex}(j,trange);
    pred=vpred{dex}(j,trange);
    pred = pred / std(xp);
    xp = xp / std(xp);
    r = corr(vtest{dex}(j,:)',vpred{dex}(j,:)');
    plot(NT*3+[1:NT],+xp-j*8, 'color', cpc(j,:));
    plot(NT*3+[1:NT],pred-j*8, 'color', cp);
    text(NT*3.6,-j*8+8,sprintf('r^2=%1.2f',r^2),'units','data','fontsize',8,...
        'fontweight','bold','color',cpc(j,:));
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
msize = 13;
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


pos(1) = h2.Position(1)+h2.Position(3)*.95;
% ------------- ball array of neurons ---------%
xb0 = pos(1)+.015;
yb0 = pos(2)-.0;
pos(3)=.08;
axes('position',[xb0 yb0 pos(3:4)]);
hold all;
ny = 10;
for j = 1
    for k = 1:ny
        xb=(j-1);
        yb=(k-1);
        if j==1
            for m = 1:nPCs
                yp = (ny-1)/3 * (m-1);
                plot([-4 xb],[yp yb],'k','linewidth',1);
            end
            %text(0.,3,'neurons','fontsize',12);
        end
        plot(xb,yb,'ko','markerfacecolor','w','linewidth',2,'color','k',...
            'markersize',7);
        
    end
end
axis tight;
axis off;

% ------------ FACE (1D) ------------%
i=i+1;
hs{i} = my_subplot(5,4,9,[.6 .6]);
hs{i}.Position(2) = hs{i}.Position(2)+.02;
hold all;
plot([0 15],[0 15],'k');
for j = 1:ndat
    plot(expv_behavior(1,j,7), expv_behavior(1,j,9),'wo',...
        'markerfacecolor',cdat(j,:),'markersize',5);
end
axis([0 .06 0 .06]);
ylabel('face (1D)');
xlabel('run+pupil+whisk (1D)');
box off;
axis square;
grid on;
grid minor;
grid minor;


% --------------------- DIMS VS VAR EXP -------------- %
i=i+1;
hs{i} = my_subplot(5,4,13,[.7 .7]);
hs{i}.Position(1) = hs{i}.Position(1) + .04;
hs{i}.Position(2) = hs{i}.Position(2)+.01;
%hs{i}.Position(2) = hs{i-1}.Position(2);
for j = 1:ndat
    semilogx(ndims0,expv_behavior(:,j,9),'color',cdat(j,:));
    hold all;
end
box off;
set(gca,'xtick',[1 4 16 64]);
xlabel('dimensions');
ylabel({'variance explained','(test set)'});
axis tight;
ylim([0 .15]);
xlim([0 128]);
hold all;
semilogx(ndims0,mean(expv_behavior(:,:,9),2),'k','linewidth',2);
axis square;
grid on;
grid minor;
grid minor;

% ------------ RASTERS
trange = [trange(5:end-13) 830+[1:70] 1200+[1:50]];
i=i+1;
hs{i}=my_subplot(5,2,6,[.8 .8]);
hs{i}.Position(1) = .34;
hs{i}.Position(3) = .63;
%hs{i}.Position(2) = hs{i-2}.Position(2)-.02;
%hs{i}.Position(1) = hs{4}.Position(1);
%hs{i}.Position(2) = hs{i}.Position(2) + .025;
%hs{i}.Position(4) = .22;
yt = ytest{dex}(isortembed{dex},trange);
ytstd = max(1e-3,std(yt,1,2));
yt = yt ./ ytstd;
imagesc(my_conv2(yt, 6, 1),[0 .5])
colormap(hs{i},flipud(cgray));
axis off;
text(0,1.15,'Neural activity (test set) sorted by 1D embedding','fontangle','normal','fontsize',8);
NN=size(yt,1);
hold all;
plot([0 30/1.2], (NN+100)*[1 1],'k');
plot([0 0]-0, NN+[0 -1000],'k');
ht=text(-.04,0,'1000 neurons','fontsize',8,'fontangle','normal');
ht.Rotation=90;
text(0,0,'    30 s','fontsize',8,'fontangle','normal');
axis tight;

hp=my_subplot(5,2,8,[.8 .8]);
hp.Position(1) = .34;
hp.Position(3) = .63;
%hp.Position(1) = hs{4}.Position(1);
%hp.Position(2) = hp.Position(2) - .00;
%hp.Position(4) = .22;
yp = ypred{dex}(isortembed{dex},trange);
yp = my_conv2(yp, 1, 1);
yp = yp ./ ytstd; % scale prediction by same amount
imagesc(yp,[0 .5])
colormap(hp,flipud(cgray));
axis off;
text(0,1.15,'Neural activity prediction (test set) from faces','fontangle','normal','fontsize',8);
ht=text(-.04,0,'1000 neurons','fontsize',8,'fontangle','normal');
ht.Rotation=90;
hold all;
plot([0 30/1.2], (NN+100)*[1 1],'k');
plot([0 0]-0, NN+[0 -1000],'k');
text(0,0,'    30 s','fontsize',8,'fontangle','normal');
axis tight;



xh=.6;
yh=.6;
%-------------- PEER PRED -----------------%
i=i+1;
hs{i}=my_subplot(5,4,17,[xh yh]);
ev = [max(expv_neurons)' (expv_behavior(6,:,9))'];
hold all;
for j = 1:size(ev,1)
    plot(ev(j,1), ev(j,2), 'wo',...
        'markerfacecolor',cdat(j,:),'markersize',5);
end
plot([0 30],[0 30],'k','linewidth',1)
axis square;
box off;
ylabel('weights');
ylim([0 .30])
xlim([0 .30]);
ylabel({'face (16D)'});
xlabel({'other neurons'});
%ht=title({'% variance explained',''},'fontsize',8);
%ht.Position(2) = ht.Position(2)-4;
grid on;
grid minor;
grid minor;


% ------------ FACE (16D) vs arousal ------------%
i=i+1;
hs{i} = my_subplot(5,4,18,[xh yh]);
hold all;
plot([0 15],[0 15],'k');
for j = 1:ndat
    plot(expv_behavior(3,j,7), expv_behavior(6,j,9),'wo',...
        'markerfacecolor',cdat(j,:),'markersize',5);
end
axis([0 .15 0 .15]);
ylabel('face (16D)');
xlabel('run+pupil+whisk (3D)');
box off;
axis square;
grid on;
grid minor;
grid minor;


% ------------ FACE (16D) + arousal ------------%
i=i+1;
hs{i} = my_subplot(5,4,19,[xh yh]);
hold all;
plot([0 15],[0 15],'k');
for j = 1:ndat
    plot(expv_behavior(6,j,10), expv_behavior(6,j,9),'wo',...
        'markerfacecolor',cdat(j,:),'markersize',5);
end
axis([0 .15 0 .15]);
ylabel('face (16D)');
xlabel('face+3D arousal (16D)');
box off;
axis square;
grid on;
grid minor;
grid minor;



% ---------- TIMELAG --------%
i=i+1;
hs{i} = my_subplot(5,4,20,[xh yh]);
tlag = tdelay*.4;
for j = 1:ndat
    plot(tlag, squeeze(expv_tlag(6,j,:)),'color',cdat(j,:));
    hold all;
end
plot(tlag,squeeze(mean(expv_tlag(6,:,:),2)),'k','linewidth',2);
axis tight;
xlabel('time from behavior (s)');
ylabel({'variance explained'});
%axis([-7*1.2 10 0 13])
ylim([0 .15]);
box off;
axis square;
grid on;
grid minor;
grid minor;



% -------------- LETTERS
hp=.065;
hy=1.22;
deffont=8;
for j = [1:length(hs)]
    if j<5
        hp0 =0.03;
        if j== 4
            hy0 = 1.3;
        else
            hy0 = 1.1;
        end
    %elseif j==5
    %    hy0=1.1;
    %    hp0=.065;
    elseif j==6
        hp0=.09;
        hy0=hy;
    elseif j==7
        hp0=.03;
        hy0=hy;
    else
        hp0=hp;
        hy0 = hy;
    end
    hpos = hs{j}.Position;
    axes('position', [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01]);
    text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
end

tl=squeeze(mean(expv_tlag(6,:,:),2));

tli = interp1(tlag,tl,[-8:.01:8]);

[~,ix1]=min(abs(max(tli)/2 - tli(1:800)));
[~,ix2]=min(abs(max(tli)/2 - tli(801:end)));

fwhm = (800-ix1+ix2)*.01

