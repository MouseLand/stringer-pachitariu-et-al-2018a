function fig4new(matroot)

load(fullfile(matroot,'faceSpectrum.mat'));
load(fullfile(matroot,'stimvar.mat'));
load(fullfile(matroot,'stimfaceRepresentation.mat'));
load(fullfile(matroot,'exResponses.mat'));

%%
default_figure([1 1 8 8.5]);

%%
rng('default');
kt = 1; % face prediction
tpts = .645e4 + .03e4 + [1:.85e4];%2.09e4;

dx = .25;
dy = .5;

sdiff = diff(stimtpt);
son = find(sdiff==1);
soff = find(sdiff==-1);

cs = colormap('cool');
cs = cs(round(linspace(36,64,5)),:);
rng(19);
cs = max(0, cs - .1*(rand(5,1)<.5));
cs = cs(randperm(size(cs,1)),:);

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
ylabel({'spont periods'},'fontsize',8);
xlabel('stim periods','fontsize',8);
%title('% variance of PCs','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

i=i+1;
hs{i} = my_subplot(4,6,5,[.55 .45]);
hs{i}.Position(2) = hs{i}.Position(2) +.02;
hs{i}.Position(1) = hs{i}.Position(1) +.01;
hold all;
for j = 1:size(Vshared,2)
    plot(Vshared(1:10,j,1)*100,'.-','markersize',8,'color',cm(j,:),'linewidth',.25)
end
box off;
ylabel({'% stim variance'},'fontsize',8);
xlabel('dimension','fontsize',8);
%title('stim-face subspace','fontweight','normal','fontangle','italic','fontsize',10)
axis square;
axis tight;
ylim([0 15]);
set(gca,'xtick',[1 5 10]);


i=i+1;
hs{i} = my_subplot(4,6,6,[.55 .45]);
hs{i}.Position(2) = hs{i}.Position(2) +.02;
hs{i}.Position(1) = hs{i}.Position(1) +.01;
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
ylabel('fraction','fontsize',8);
xlabel('neural weights','fontsize',8);
%title('1st PC weights','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

i=i+1;
hs{i} = my_subplot(4,2,4,[1.2 1.2]);
%hs{i}.Position(1)=hs{i-3}.Position(1);
hs{i}.Position(2)=hs{i}.Position(2)+.05;
hs{i}.Position(1)=hs{i}.Position(1)+.03;
axis off;
im=imread('schematicNEW.png');
image(im);
axis tight;
axis image
axis off;

kt = 1; % use face projections
tshared1 = tshared{kt} * sign(mean(Ushared{1,kt}));
tshared2 = tshared{2} * sign(mean(Ushared{2,2}));
i=i+1;
hs{i} = my_subplot(2,2,3,[.85 .75]);
hs{i}.Position(1) = hs{1}.Position(1);
hs{i}.Position(2) = hs{i}.Position(2)+.03;
%hs{i}.Position(3) = .5;
dobin=1;
tproj1=tprojS{1}(1:3,:);
tproj1= sign(skewness(tproj1,1,2)) .* tproj1;
tproj2=tprojS{2}(1:3,:);
tproj2= sign(skewness(tproj2,1,2)) .* tproj2;
titles = {'stim-only','spont-only','face-only'};
fig_traces(sprojS{kt}(1:2:10,tpts),tshared1(tpts),tproj1(:,tpts),tshared2(tpts),tproj2(:,tpts),cs,dobin,titles);
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


dx = .7;
dy = .65;
%istims1 = istims;
%istims2 = istims;
%istims1(1:end/2) = 0;
%istims2(floor(length(istims2)/2)+1:end) = 0;

i=i+1;
hs{i}=my_subplot(6,4,15,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.005;
hold all;
k=0;
isti = [22 1];
hold all;
k=0;
cs(1,:) = [.5 .2 .7];
cs(2,:) = [1 .5 1];
for j = isti
    k=k+1;
    plot(projstim{1}{3}(istims{1}==j,isti(1)),projstim{1}{3}(istims{1}==j,isti(2)),'.','color',cs(k,:));
	plot(Rfit{1}{3}(istims{1}==j,isti(1)), Rfit{1}{3}(istims{1}==j,isti(2)), 'k','linewidth',0.5);
    %text(1.1,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
end
text(1.2,1,'stim 1', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.2,0.85,'stim 2', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.2,0.7,{'multiplicative','model'},'HorizontalAlignment','right','fontangle','normal','fontsize',8);
axis tight;
axis([0 250 0 220]);
axis square;
box off;
xlabel('stim-only 1');
ylabel('stim-only 2');

i=i+1;
hs{i}=my_subplot(6,4,16,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.005;
id = [3 1 2];
k=1;
normv = zeros(4,4);
normv0= sum(vsigstimspont(:,:,id(k),:),1)./sum(vsigstimspont(:,2,id(k),:),1);
normv0 = squeeze(normv0);
normv([1 4],:) = normv0([1 3],:);
normv(2:3,:) = squeeze(fitmult(:,id(k),:));

for d = 1:4
	hold all;
	for j = 1:3
		plot([1:4]+(k-1)*3, normv(:,d),'o','color',cm(d,:),'markersize',4);
	end
end
title('stim-only dims','fontweight','normal');
set(gca,'xtick',[1:4],'xticklabel',{'avg. model','mult. model','affine model','spont period var','signal var','stim period','spont period','signal var','stim period','spont period'});
set(gca,'XTickLabelRotation',45);
xlim([.5 4.5]);
ylabel({'fraction of variance'});
axis square


% normalize projstim by dimensions
i=i+1;
hs{i}=my_subplot(6,2,10,[.8 .4]);
hs{i}.Position(2)=hs{i}.Position(2)-.01;
x = multgain{1}{3};
y = projstim{1}{4};
a = (x'*x)\(x'*y);
hold all;
tstim=[100:500];
plot(y(tstim),'color',.4*[1 1 1],'linewidth',1);
plot(x(tstim)*a,'r','linewidth',0.5)
text(0,1.25,'stim-face shared dim','color',.4*[1 1 1],'fontangle','normal','fontsize',8);
text(0,1,'multiplicative gain','color','r','fontangle','normal','fontsize',8);
box off;
axis tight;
axis off;

i=i+1;
hs{i}=my_subplot(6,4,19+4,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.05;
hold all;
k=0;
for j = isti
    k=k+1;
    plot(projstim{1}{1}(istims{1}==j,1),projstim{1}{1}(istims{1}==j,2),'.','color',cs(k,:));
	%plot(Rfit{1}{1}(istims{1}==j,isti(1)), Rfit{1}{1}(istims{1}==j,isti(2)), 'k','linewidth',0.5);
    %text(1.3,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
end
axis tight;
axis square;
box off;
xlabel('face-only 1');
ylabel('face-only 2');

i=i+1;
hs{i}=my_subplot(6,4,24,[dx dy]);
hs{i}.Position(2)=hs{i}.Position(2)+0.05;
id = [1];
k=1;
normv = zeros(4,4);
normv0= sum(vsigstimspont(:,:,id(k),:),1)./sum(vsigstimspont(:,2,id(k),:),1);
normv0 = squeeze(normv0);
normv([1 4],:) = normv0([1 3],:);
normv(2:3,:) = squeeze(fitmult(:,id(k),:));

for d = 1:4
	hold all;
	for j = 1:3
		plot([1:4]+(k-1)*3, normv(:,d),'o','color',cm(d,:),'markersize',4);
	end
end
title('face-only dims','fontweight','normal');
set(gca,'xtick',[1:4],'xticklabel',{'avg. model','mult. model','affine model','spont period var','signal var','stim period','spont period','signal var','stim period','spont period'});
set(gca,'XTickLabelRotation',45);
ylim([0 1.05]);
xlim([.5 4.5]);
ylabel({'fraction of variance'});
axis square;

tstr{1} = 'Behavior and neural activity with/without visual stimulation';
tstr{2} = 'Variance of face PCs';
tstr{3} = 'Stim-face shared dims';
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
    if j ==1 || j==6
        hp0 =.04;
        hy0 = 1.06;
    elseif j==5 
        hy0=.95;
        hp0 = -0.08;
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


%% plots of variance as a function of dimension

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
