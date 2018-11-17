function fig4new(matroot)

load(fullfile(matroot,'faceSpectrum.mat'));
load(fullfile(matroot,'stimvar.mat'));
load(fullfile(matroot,'stimfaceRepresentation.mat'));
load(fullfile(matroot,'exResponses.mat'));

%%
default_figure([1 1 8 8.5]);

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
hs{i} = my_subplot(2,2,1,[.85 .8]);
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
ht=text(-.08,.75,'Face motion energy PCs','fontsize',8,'fontangle','normal','color',cpc(5,:),...
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
ht=text(-.08,.25,'Example neurons','fontsize',8,'fontangle','normal',...
    'HorizontalAlignment','center');
ht.Rotation=90;

axes('position',[hs{i}.Position(1) hs{i}.Position(2)-yh hs{i}.Position(3) yh*.8]);
fig_stimon(stimtpt(tpts))

cm = colormap('hot');
cm = cm(32+[1:4:4*4],:);

i=i+1;
hs{i} = my_subplot(4,6,4,[.55 .45]);
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
ylabel({'spont periods'},'fontsize',8);
xlabel('stim periods','fontsize',8);
%title('% variance of PCs','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

i=i+1;
hs{i} = my_subplot(4,6,5,[.55 .45]);
hs{i}.Position(2) = hs{i}.Position(2) +.02;
hs{i}.Position(1) = hs{i}.Position(1) +.01;
hold all;
for d = 1:size(Vshared,2)
    plot(Vshared(1:10,d,ktype)*100,'.-','markersize',8,'color',cm(d,:),'linewidth',.25)
end
box off;
ylabel({'% stim variance'},'fontsize',8);
xlabel('dimension','fontsize',8);
%title('stim-face subspace','fontweight','normal','fontangle','italic','fontsize',10)
axis square;
axis tight;
ylim([0 14]);
set(gca,'xtick',[1 5 10]);


i=i+1;
hs{i} = my_subplot(4,6,6,[.55 .45]);
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
ylabel('fraction','fontsize',8);
xlabel('neural weights','fontsize',8);
%title('1st PC weights','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

i=i+1;
hs{i} = my_subplot(4,2,4,[1.2 1.2]);
%hs{i}.Position(1)=hs{i-3}.Position(1);
hs{i}.Position(2)=hs{i}.Position(2)+.065;
hs{i}.Position(1)=hs{i}.Position(1)+.03;
axis off;
im=imread('schematicNEW.png');
image(im);
axis tight;
axis image
axis off;

tshared1 = tshared{1} * sign(mean(Ushared{1,1}));
tshared2 = tshared{2} * sign(mean(Ushared{2,2}));
i=i+1;
hs{i} = my_subplot(2,2,3,[.85 .82]);
hs{i}.Position(1) = hs{1}.Position(1);
hs{i}.Position(2) = hs{i}.Position(2)+.05;
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
text(.2,1.0,'running','color',[.2 .8 .2],'fontangle','normal','fontsize',8);
hold all;
plot(1 + [0 60*5],-6 * [1 1],'k','linewidth',2);
text(0.06,-.1,'5 min','horizontalalignment','center','fontsize',8,'fontangle','normal');
axis off;
axis tight;


dx = .7;
dy = .65;
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
text(xy,1.15-xx,'stim-only','color',cp(3,:),'fontangle','normal','fontsize',8);
text(xy,1-xx,'behav-only','color',cp(1,:),'fontangle','normal','fontsize',8);
text(xy,.85-xx,'spont-only','color',cp(2,:),'fontangle','normal','fontsize',8);
text(xy,.7-xx,'stim-behav','color',cp(4,:),'fontangle','normal','fontsize',8);
xlabel('spont period');
ylabel('stim period');
title('projection variances','fontweight','normal');
set(gca,'ytick',10.^[2:4]);
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
    %text(1.1,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
end
text(1.2,1,'stim 1', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.2,0.85,'stim 2', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.2,0.7,{'multiplicative','gain model'},'color','r','HorizontalAlignment','right','fontangle','normal','fontsize',8);
%axis tight;
axis([-50 180 -20 150]);
%axis([0 220 0 250]);
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
    %text(1.3,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
end
axis tight;
axis square;
box off;
xlabel('behav-only proj1');
ylabel('behav-only proj2');
text(1.4,1,'stim 1', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.4,0.85,'stim 2', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);


% i=i+1;
% hs{i}=my_subplot(6,4,20,[dx dy]);
% hs{i}.Position(2)=hs{i}.Position(2)+0.02;
% k=1;
% normv = squeeze(vsigstimspont(4,[3 2 1 4],:));% ./ sum(vsigstimspont(2,[1 3:4],:),2));
% randv = mean(squeeze(mean(vsigstimspont(4,5:end,:))));
% for d = 1:4
% 	hold all;
% 	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),'o-','color',cm(d,:),'markersize',4);
% end
% plot([0 size(normv,1)+1], randv * [1 1],'k--','linewidth',2);
% text(1.3,.22,'random','HorizontalAlignment','right','fontsize',8)
% %title({'stim responses','of projections'},'fontweight','normal');
% set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'stim-only','spont-only','behav-only','stim-behav'});
% set(gca,'XTickLabelRotation',30);
% xlim([.5 size(normv,1)+.5]);
% ylabel({'signal-to-noise ratio'});
% axis square;


i=i+1;
hs{i}=my_subplot(6,4,20,[dx dy]);
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
ylim([.0 1.0]);
xlim([.5 size(normv,1) + .5]);
ylabel({'variance explained'});
axis square;

% mult gain and face pred
i=i+1;
d=1;
hs{i}=my_subplot(6,2,12,[.78 .62]);
hs{i}.Position(2)=hs{i}.Position(2)-.01;
hs{i}.Position(1)=hs{i}.Position(1)-.02;
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
tstim=300:800;
r = corr(y(ttest),mg(ttest));
plot(mg(ttest(tstim)),'r','linewidth',1)
plot(y(ttest(tstim))*as,'color',.6*[1 1 1],'linewidth',1);
plot(fg(:,ttest(tstim))'*a,'b','linewidth',0.5)
text(0,1.15,'gain','color','r','fontangle','normal','fontsize',8);
text(.25,1.15,sprintf('face prediction of gain, r=%1.2f',facegain_r2(1)),'color','b','fontangle','normal','fontsize',8);
text(.3,0.1,sprintf('stim-behav shared dim, r=%1.2f',r),'color',.4*[1 1 1],'fontangle','normal','fontsize',8);
box off;
axis tight;
axis off;
plot([0 60/1.2],[1 1]*-110,'k','linewidth',2)
text(0,-.02,'1 min','fontangle','normal','fontsize',8);


%normv(2:3,:) = squeeze(fitmult(:,id(k),:)) .* normv(4,:);


tstr{1} = 'Behavior and neural activity with/without visual stimulation';
tstr{2} = 'Variance of face PCs';
tstr{3} = {'Stim-behav shared','         subspace'};
tstr{4} = 'Top dimension';
tstr{5} = '';
tstr{6} = 'Projections of neural activity';
for k = 7:length(hs)
    tstr{k} = '';
end

% -------------- LETTERS
hp=.05;
hy=1.22;
deffont=8;
for j = [1:length(hs)]
    if j ==1
        hp0 =.04;
        hy0 = 1.06;
	elseif j==6
		hp0 =.04;
        hy0 = 1.03;
    elseif j==5 
        hy0=.95;
        hp0 = -0.08;
    elseif j==9
		hy0=1.2;
		hp0=.03;
	elseif j>6
        hy0=hy;
        hp0=.03;
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
%%
print(fullfile(matroot, 'fig4new.pdf'),'-dpdf');


%% plots of variance as a function of dimension
clf;
cm=[0 0 1; 0 .5 .5; 1 0 .5; .5 .5 .5];
for k=1:4
	plot(squeeze(vsigstimspont(3,k,:)),squeeze(vsigstimspont(2,k,:)),'o','color',cm(k,:))
	hold all;
end
hold all;
plot([25 5e3],[25 5e3],'k')
xlabel('spont periods');
ylabel('stim periods');
axis tight;

%%
% dx = .7;
% dy = .65;
% i=i+1;
% cm = [.4 .4 .4; .8 .5 .8; .4 .4 1];
% hs{i}=my_subplot(6,4,15,[dx dy]);
% hs{i}.Position(2)=hs{i}.Position(2)+.005;
% hold all;
% for k = 1:3
%     plot(vsigstimspont(:,k,3,1),'color',cm(k,:))
% end
% axis square;
% axis tight;
% box off;
% ylabel('variance');
% xlabel('dimension');
% ylim([0 max(reshape(vsigstimspont(:,:,3,1),[],1))])
% text(1.2,.6,'signal variance','color',cm(1,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
% text(1.2,.9,'stim periods','color',cm(2,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
% text(1.2,.75,'spont periods','color',cm(3,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
% text(0.5,1.1,'stim-only','fontweight','normal','HorizontalAlignment','center','FontAngle','normal');
% 
% i=i+1;
% hs{i}=my_subplot(6,4,15+4,[dx dy]);
% hs{i}.Position(2)=hs{i}.Position(2)+.005;
% hold all;
% for k = 1:3
%     plot(vsigstimspont(:,k,1,1),'color',cm(k,:))
% end
% axis square;
% axis tight;
% box off;
% text(0.5,1.1,'face-only','fontweight','normal','HorizontalAlignment','center','FontAngle','normal');
% 
% i=i+1;
% hs{i}=my_subplot(6,4,15+8,[dx dy]);
% hs{i}.Position(2)=hs{i}.Position(2)+0.005;
% hold all;
% for k = 1:3
%     plot(vsigstimspont(:,k,2,1),'color',cm(k,:))
% end
% axis square;
% axis tight;
% box off;
% text(0.5,1.1,'spont-only','fontweight','normal','HorizontalAlignment','center','FontAngle','normal');


% 
% i=i+1;
% cm = [1 .4 .5; .6 0 .6; .3 .3 1];
% hs{i}=my_subplot(6,4,16,[dx dy]);
% hs{i}.Position(2)=hs{i}.Position(2)+.005;
% hold all;
% for k = 1:3
%     shadedErrorBar([1:32],mean(vsigstimspont(:,k,3,:),4),std(vsigstimspont(:,k,3,:),1,4)/sqrt(3),{'color',cm(k,:)})
% end
% text(1.2,.6,'signal variance','color',cm(1,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
% text(1.2,.9,'stim periods','color',cm(2,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
% text(1.2,.75,'spont periods','color',cm(3,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
% axis square;
% axis tight;
% box off;
% ylabel('variance');
% xlabel('dimension');
% text(0.5,1.1,'stim-only','fontweight','normal','HorizontalAlignment','center','FontAngle','normal');
% ylim([0 max(reshape(vsigstimspont(:,:,3,1),[],1))])
% 
% i=i+1;
% hs{i}=my_subplot(6,4,16+4,[dx dy]);
% hs{i}.Position(2)=hs{i}.Position(2)+.005;
% hold all;
% for k = 1:3
%     shadedErrorBar([1:32],mean(vsigstimspont(:,k,1,:),4),std(vsigstimspont(:,k,1,:),1,4)/sqrt(3),{'color',cm(k,:)})
% end
% axis square;
% axis tight;
% box off;
% text(0.5,1.1,'face-only','fontweight','normal','HorizontalAlignment','center','FontAngle','normal');
% 
% i=i+1;
% hs{i}=my_subplot(6,4,16+8,[dx dy]);
% hs{i}.Position(2)=hs{i}.Position(2)+.005;
% hold all;
% for k = 1:3
%     shadedErrorBar([1:32],mean(vsigstimspont(:,k,2,:),4),std(vsigstimspont(:,k,2,:),1,4)/sqrt(3),{'color',cm(k,:)})
% end
% axis square;
% axis tight;
% box off;
% text(0.5,1.1,'spont-only','fontweight','normal','HorizontalAlignment','center','FontAngle','normal');
