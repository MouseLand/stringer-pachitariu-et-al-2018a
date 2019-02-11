function fig4new(matroot)

load(fullfile(matroot,'stimvar_with_decode.mat'));
load(fullfile(matroot,'faceSpectrum.mat'));
load(fullfile(matroot,'stimfaceRepresentation.mat'));
load(fullfile(matroot,'exResponses.mat'));

%%
close all;
default_figure([1 1 7.25 7.25]);

%%
rng('default');
ktype = 1; % spont prediction
tpts = .645e4 + .03e4 + [1:.85e4];%2.09e4;

dx = .25;
dy = .5;

sdiff = diff(stimtpt);
son = find(sdiff==1);
soff = find(sdiff==-1);

cs = [];
cs(1,:) = [.5 .2 .7];
cs(2,:) = [1 .3 1];
cs(3,:) = [.8 .3 .8];
cs(4,:) = [.5 0 .5];
cs(5,:) = [.7 .4 .7];

yh=.02;
X = [.05 .2 .35 .72];
Y = [.05 .05+dy*.55 .64];

clf;

clear hs;
i=0;

i=i+1;
hs{i} = my_subplot(2,2,1,[.9 .8]);
hs{i}.Position(1) = hs{i}.Position(1) -.0;
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
ht=text(-.04,.75,'Face motion energy PCs','fontsize',6,'fontangle','normal','color',cpc(5,:),...
    'HorizontalAlignment','center');
ht.Rotation=90;
rng(12)
%ineur = [100 20 1999 1800 1900 1902 1998 1990 1997 1996];
ineur = [1:1000:10000];
for j = 1:10 %size(Fex,2)
    fplot = Fex(ineur(j),tpts);
    fplot = bin2d(fplot(:),4);
    plot(tpts(1:4:end),zscore(fplot)/4-j*2.5,'color','k','linewidth',.5)%cn(j,:))
end
axis off;
axis tight;
plot((find(diff(stimtpt(tpts))==1,1)+tpts(1)-1) * [1 1], [-10*2.5 10.5*2.5],'k--','color',.5*[1 1 1]);
plot((find(diff(stimtpt(tpts))==-1,1)+tpts(1)-1) * [1 1], [-10*2.5 10.5*2.5],'k--','color',.5*[1 1 1]);
ht=text(-.04,.25,'Example neurons','fontsize',6,'fontangle','normal',...
    'HorizontalAlignment','center');
ht.Rotation=90;
set(gca,'fontsize',6);

axes('position',[hs{i}.Position(1) hs{i}.Position(2)-yh hs{i}.Position(3) yh*.8]);
fig_stimon(stimtpt(tpts))
set(gca,'fontsize',6);

cm = colormap('hot');
cm = cm(32+[1:4:4*4],:);

xh=.5;
yh=.45;

i=i+1;
hs{i} = my_subplot(4,6,4,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) +.02;
loglog([1e-3 1],[1e-3 1],'k','linewidth',.5)
hold all;
for d = 1:4
    x = vface(1:10,1,d)/nansum(vface(:,1,d))*1;
    y = vface(1:10,2,d)/nansum(vface(:,2,d))*1;
    loglog(x,y,'.','color',cm(d,:),'markersize',8)
end
%ylim([.8 1.2]);
box off;
axis tight;
set(gca,'ytick',10.^[-2:0],'yticklabel',{'0.01','0.1','1'},...
    'xtick',10.^[-2:0],'xticklabel',{'0.01','0.1','1'});
ylabel({'spont periods'},'fontsize',6);
xlabel('stim periods','fontsize',6);
%title('% variance of PCs','fontweight','normal','fontangle','italic','fontsize',10)
axis square;
set(gca,'fontsize',6);

i=i+1;
hs{i} = my_subplot(4,6,5,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) +.02;
hs{i}.Position(1) = hs{i}.Position(1) +.01;
hold all;
for d = 1:size(Vshared,2)
    plot(Vshared(1:10,d,ktype)*100,'.-','markersize',8,'color',cm(d,:),'linewidth',.25)
end
box off;
ylabel({'% stim variance'},'fontsize',6);
xlabel('dimension','fontsize',6);
%title('stim-face subspace','fontweight','normal','fontangle','italic','fontsize',10)
axis square;
axis tight;
ylim([0 14]);
set(gca,'xtick',[1 5 10]);
set(gca,'fontsize',6);

i=i+1;
hs{i} = my_subplot(4,6,6,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) +.02;
hs{i}.Position(1) = hs{i}.Position(1) +.01;
hold all;
for d = 1:size(Ushared,1)
    uplot = Ushared{d,ktype};
    uplot = uplot * sign(mean(uplot));
    bed = [-.02:.001:.035];
    histogram(uplot,20,'binedges',bed,'edgecolor',cm(d,:),...
        'normalization','probability','displaystyle','stairs');
	upos(d) = mean(Ushared{d,ktype}>0);
end
box off;
axis tight;
xlim([-.01 .03]);
ylabel('fraction','fontsize',6);
xlabel('neural weights','fontsize',6);
%title('1st PC weights','fontweight','normal','fontangle','italic','fontsize',10)
axis square;
set(gca,'fontsize',6);

i=i+1;
hs{i} = my_subplot(4,4,7,[1.1 1.1]);
%hs{i}.Position(1)=hs{i-3}.Position(1);
hs{i}.Position(2)=hs{i}.Position(2)+.065;
hs{i}.Position(1)=hs{i}.Position(1)+.03;
axis off;
im=imread('schematicNEW.png');
image(im);
axis tight;
axis image
axis off;
set(gca,'fontsize',6);

yh=.02;

dx = .7;
dy = .65;

i=i+1;
hs{i}=my_subplot(6,4,8,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)-0.06;
hs{i}.Position(1)=hs{i}.Position(1)+0.03;
k=1;
normv = squeeze(decoding([3 2 1],:));% ./ sum(vsigstimspont(2,[1 3:4],:),2));
randv = mean(decoding(5:end,:),1);
normv = 1 - [normv(1:3, :); randv];% normv(end,:)];
for d = 1:4
	hold all;
	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),'o-','color',cm(d,:),'markersize',4);
end
%title({'stim responses','of projections'},'fontweight','normal');
set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'stim-only','spont-only','behav-only','random 32D', 'stim-behav'});
set(gca,'XTickLabelRotation',30);
xlim([.5 size(normv,1)+.5]);
ylabel({'decoding error'});
axis square;
set(gca,'fontsize',6);


tshared1 = tshared{1} * sign(mean(Ushared{1,1}));
tshared2 = tshared{2} * sign(mean(Ushared{2,2}));
i=i+1;
hs{i} = my_subplot(2,2,3,[.9 .87]);
hs{i}.Position(1) = hs{1}.Position(1);
hs{i}.Position(2) = hs{i}.Position(2)+.055;
%hs{i}.Position(3) = .5;
dobin=1;
tproj1=tprojS{1}(1:3,:);
tproj1= sign(skewness(tproj1,1,2)) .* tproj1;
tproj2=tprojS{2}(1:3,:);
tproj2= sign(skewness(tproj2,1,2)) .* tproj2;
titles = {'stim-only','behav-only','spont-only'};
fig_traces(sprojS{ktype}(1:2:10,tpts),tshared1(tpts),tproj1(:,tpts),tshared2(tpts),tproj2(:,tpts),cs,dobin,titles);
axes('position',[hs{i}.Position(1) hs{i}.Position(2)-yh hs{i}.Position(3) yh*.8]);
fig_stimon(stimtpt(tpts),0);
axes('position',[hs{i}.Position(1) hs{i}.Position(2)-4*yh hs{i}.Position(3) yh*2.8]);
plot(bin2d(runspeed(tpts),4),'color',[.2 .8 .2]);
text(.2,1.0,'running','color',[.2 .8 .2],'fontangle','normal','fontsize',6);
hold all;
plot(1 + [0 60*5/1.3],-6 * [1 1],'k','linewidth',2);
text(0.06,-.03,'5 min','horizontalalignment','center','fontsize',6,'fontangle','normal');
axis off;
axis tight;
set(gca,'fontsize',6);
%istims1 = istims;
%istims2 = istims;
%istims1(1:end/2) = 0;
%istims2(floor(length(istims2)/2)+1:end) = 0;



i=i+1;
hs{i}=my_subplot(6,4,15,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.025;
cp=[.4 .6 1; 0 0 0; 0.8 0.3 0.8; .5 .5 .5];
ll=[squeeze(vsigstimspont(3,:,:)),squeeze(vsigstimspont(2,:,:))];
ll = ll(:);
loglog([min(ll) max(ll)],[min(ll) max(ll)],'k','linewidth',.5)
hold all;
for k=1:4
	loglog(squeeze(vsigstimspont(3,k,:)),squeeze(vsigstimspont(2,k,:)),'x','color',cp(k,:))
	hold all;
end
xx = .52;
xy = .65;
text(xy,1.15-xx,'stim-only','color',cp(3,:),'fontangle','normal','fontsize',6);
text(xy,1-xx,'behav-only','color',cp(1,:),'fontangle','normal','fontsize',6);
text(xy,.85-xx,'spont-only','color',cp(2,:),'fontangle','normal','fontsize',6);
text(xy,.7-xx,'stim-behav','color',cp(4,:),'fontangle','normal','fontsize',6);
xlabel('spont period');
ylabel('stim period');
title('projection variances','fontweight','normal');
set(gca,'ytick',10.^[2:4]);
set(gca,'fontsize',6);
axis square;
box off;
axis tight;


isti = [22 1];
cs(1,:) = [.5 .2 .7];
cs(2,:) = [1 .3 1];

i=i+1;
hs{i}=my_subplot(6,4,16,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.025;
hold all;
k=0;
for j = isti
    k=k+1;
    plot(projstim{1}{3}(istims{1}==j,isti(1)),projstim{1}{3}(istims{1}==j,isti(2)),...
		'.','color',cs(k,:),'markersize',5);
	plot(Rfit{1}{3}(istims{1}==j,isti(1)), Rfit{1}{3}(istims{1}==j,isti(2)), 'r','linewidth',0.5);
    %text(1.1,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',6,'HorizontalAlignment','right')
end
text(1.2,1,'stim 1', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',6);
text(1.2,0.85,'stim 2', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',6);
text(1.2,0.7,{'multiplicative','gain model'},'color','r','HorizontalAlignment','right','fontangle','normal','fontsize',6);
%axis tight;
axis([-50 180 -20 150]);
%axis([0 220 0 250]);
set(gca,'fontsize',6);
axis square;
box off;
xlabel('stim-only proj1');
ylabel('stim-only proj2');

i=i+1;
hs{i}=my_subplot(6,4,19,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.02;
hold all;
k=0;
for j = isti
    k=k+1;
    plot(projstim{1}{ktype}(istims{1}==j,1),projstim{1}{ktype}(istims{1}==j,2)*-1,...
		'.','color',cs(k,:),'markersize',5);
	%plot(Rfit{1}{1}(istims{1}==j,isti(1)), Rfit{1}{1}(istims{1}==j,isti(2)), 'k','linewidth',0.5);
    %text(1.3,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',6,'HorizontalAlignment','right')
end
axis tight;
axis square;
box off;
set(gca,'fontsize',6);
xlabel('behav-only proj1');
ylabel('behav-only proj2');
text(1.4,1,'stim 1', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',6);
text(1.4,0.85,'stim 2', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',6);


i=i+1;
hs{i}=my_subplot(6,4,20,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.02;
k=1;
normv = squeeze(vsigstimspont(4,[3 2 1 4],:));% ./ sum(vsigstimspont(2,[1 3:4],:),2));
randv = mean(squeeze(mean(vsigstimspont(4,5:end,:))));
for d = 1:4
	hold all;
	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),'o-','color',cm(d,:),'markersize',4);
end
plot([0 size(normv,1)+1], randv * [1 1],'k--','linewidth',2);
text(1.3,.22,'random','HorizontalAlignment','right','fontsize',6)
%title({'stim responses','of projections'},'fontweight','normal');
set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'stim-only','spont-only','behav-only','stim-behav'});
set(gca,'XTickLabelRotation',30);
set(gca,'fontsize',6);
xlim([.5 size(normv,1)+.5]);
ylabel({'signal-to-noise ratio'});
axis square;


i=i+1;
hs{i}=my_subplot(6,4,23,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.02;
id = [1];
normv = squeeze(vsigstimspont(1,[3],:) ./ vsigstimspont(2,[3],:));
normv = [normv'; squeeze(fitmult(:,3,:))];
for d = 1:4
	hold all;
	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),'o-','color',cm(d,:),'markersize',4);
end
%title('spont periods','fontweight','normal');
set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'avg. model','mult. model','affine model'});
set(gca,'XTickLabelRotation',30);
set(gca,'fontsize',6);
ylim([.0 1.0]);
xlim([.5 size(normv,1) + .5]);
ylabel({'variance explained'});
axis square;

% mult gain and face pred
i=i+1;
d=1;
hs{i}=my_subplot(6,4,24,[.9 .62]);
hs{i}.Position(2)=hs{i}.Position(2)-.015;
hs{i}.Position(1)=hs{i}.Position(1)-.01;
[ttrain,ttest]=splitInterleaved(length(multgain{d}{3}),30,0.5,1);
mg = multgain{d}{3};
fg = facest{d};
mg = mg-mean(mg);
fg = fg-mean(fg,2);
a=(fg(:,ttrain)*fg(:,ttrain)' + 100*eye(500))\(fg(:,ttrain)*mg(ttrain));
% x = x - mean(x);
y = projstim{1}{4};
y = y - mean(y);
% z = facest{1}';
% z = z - mean(z,1);
as = (y(ttrain)'*y(ttrain))\(y(ttrain)'*mg(ttrain));
% a2 = (z'*z + 300*eye(size(z,2)))\(z'*x);
hold all;
ttest = find(ttest);
tstim=300:700;
r = corr(y(ttest),mg(ttest));
plot(mg(ttest(tstim)),'r','linewidth',1)
%plot(y(ttest(tstim))*as,'color',.6*[1 1 1],'linewidth',1);
plot(fg(:,ttest(tstim))'*a,'b','linewidth',0.5)
text(0,1.3,'gain','color','r','fontangle','normal','fontsize',6);
text(.25,1.3,sprintf('face prediction of gain\n          r=%1.2f',face_pred_gain_r2(1)),...
	'color','b','fontangle','normal','fontsize',6);
%text(.3,0.1,sprintf('stim-behav shared dim, r=%1.2f',r),'color',.4*[1 1 1],'fontangle','normal','fontsize',6);
box off;
axis tight;
axis off;
plot([0 60/1.2],[1 1]*-110,'k','linewidth',2)
text(0,-.02,'1 min','fontangle','normal','fontsize',6);


%normv(2:3,:) = squeeze(fitmult(:,id(k),:)) .* normv(4,:);


tstr{1} = 'Behavior and neural activity with/without visual stimulation';
tstr{2} = 'Variance of face PCs';
tstr{3} = {'Stim-behav shared','         subspace'};
tstr{4} = 'Top dimension';
tstr{5} = '';
tstr{6} = '';
tstr{7} = 'Projections of neural activity';
for k = 8:length(hs)
    tstr{k} = '';
end

% -------------- LETTERS
hp=.05;
hy=1.22;
deffont=6;
for j = [1:length(hs)]
    if j ==1
        hp0 =.02;
        hy0 = 1.06;
	elseif j==7
		hp0 =.02;
        hy0 = 1.0;
    elseif j==5 
        hy0=.85;
        hp0 = 0.005;
    elseif j==length(hs)
		hy0=1.4;
		hp0=.03;
	elseif j>7 || j==6
        hy0=hy;
        hp0=.03;
	else
        hp0=hp;
        hy0 = hy;
    end
    hpos = hs{j}.Position;
    lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
    axes('position', lpos);
    text(0,0, char(64+j),'fontsize',deffont+4,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
    
    axes('position', [lpos(1)+.02 lpos(2)-.005 lpos(3) lpos(4)]);
    text(0,0, tstr{j},'fontangle','normal','fontsize',deffont);
    axis([0 1 0 1]);
    axis off;
end
%%
print(fullfile(matroot, 'fig4new.pdf'),'-dpdf');

