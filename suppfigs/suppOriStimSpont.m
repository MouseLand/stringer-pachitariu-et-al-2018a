function fig4new(matroot)

load(fullfile(matroot,'stimvar_ori32.mat'));
load(fullfile(matroot,'faceSpectrum_ori32.mat'));

%%
close all;
default_figure([8 1 4 8]);

%%
rng('default');
ktype = 1; % spont prediction

clear hs;
i=0;
clf;
cm = colormap('hot');
cm = cm(32+[1:4:4*4],:);

dx = .5;
dy = .5;

i=i+1;
hs{i} = my_subplot(5,3,1,[dx dy]);
loglog([1e-3 1],[1e-3 1],'k','linewidth',.5)
hold all;
for d = 1:4
    x = vface(1:10,1,d)/nansum(vface(:,1,d))*1;
    y = vface(1:10,2,d)/nansum(vface(:,2,d))*1;
    loglog(x,y,'.','color',cm(d,:),'markersize',8)
end
ylim([.8 1.2]);
box off;
axis tight;
set(gca,'ytick',10.^[-2:0],'yticklabel',{'0.01','0.1','1'},...
    'xtick',10.^[-2:0],'xticklabel',{'0.01','0.1','1'});
ylabel({'spont periods'},'fontsize',8);
xlabel('stim periods','fontsize',8);
%title('% variance of PCs','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

i=i+1;
hs{i} = my_subplot(5,3,2,[dx dy]);
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
ylim([0 8]);
set(gca,'xtick',[1 5 10]);


i=i+1;
hs{i} = my_subplot(5,3,3,[dx dy]);
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
ylim([0 0.2]);
ylabel('fraction','fontsize',8);
xlabel('neural weights','fontsize',8);
%title('1st PC weights','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

dx = .55;
dy = .55;


i=i+1;
hs{i}=my_subplot(5,2,3,[dx dy]);
k=1;
normv = squeeze(decoding([3 2 1],:,2));% ./ sum(vsigstimspont(2,[1 3:4],:),2));
randv = mean(decoding(5:end,:,2),1);
normv = [normv(1:3, :); randv];% normv(end,:)];
normv = normv * 360/33;
for d = 1:4
	hold all;
	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),'o-','color',cm(d,:),'markersize',4);
end
%title({'stim responses','of projections'},'fontweight','normal');
set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'stim-only','spont-only','behav-only','random 32D'});
set(gca,'XTickLabelRotation',30);
xlim([.5 size(normv,1)+.5]);
ylabel({'error (degrees)'});
axis square;

i=i+1;
hs{i}=my_subplot(5,2,4,[dx dy]);
k=1;
normv = squeeze(decoding([3 2 1],:,1));% ./ sum(vsigstimspont(2,[1 3:4],:),2));
randv = mean(decoding(5:end,:,1),1);
normv = [normv(1:3, :); randv];% normv(end,:)];
normv = normv * 360/33;
for d = 1:4
	hold all;
	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),'o-','color',cm(d,:),'markersize',4);
end
%title({'stim responses','of projections'},'fontweight','normal');
set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'stim-only','spont-only','behav-only','random 32D'},'ytick',[0:30:90]);
set(gca,'XTickLabelRotation',30);
xlim([.5 size(normv,1)+.5]);
ylabel({'error (degrees)'});
ylim([0 90]);
axis square;

i=i+1;
hs{i}=my_subplot(5,2,5,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.025;
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


isti = [5 2];
cs(1,:) = [.5 .2 .7];
cs(2,:) = [1 .3 1];

i=i+1;
hs{i}=my_subplot(5,2,6,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.025;
hold all;
k=0;
for j = isti
    k=k+1;
    plot(projstim{1}{3}(istims{1}==j,isti(1)),projstim{1}{3}(istims{1}==j,isti(2)),...
		'.','color',cs(k,:),'markersize',5);
	plot(Rfit{1}{3}(istims{1}==j,isti(1)), Rfit{1}{3}(istims{1}==j,isti(2)), 'r','linewidth',0.5);
    %text(1.1,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
end
text(1.2,1,'11^{\circ} grating', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.2,0.8,'55^{\circ} grating', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.2,0.6,{'multiplicative','gain model'},'color','r','HorizontalAlignment','right','fontangle','normal','fontsize',8);
%axis tight;
axis([-20 180 -20 200]);
%axis([0 220 0 250]);
axis square;
box off;
xlabel('stim-only proj1');
ylabel('stim-only proj2');

i=i+1;
hs{i}=my_subplot(5,2,7,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.02;
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
text(1.6,1,'11^{\circ} grating', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.6,0.85,'55^{\circ} grating', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);


i=i+1;
hs{i}=my_subplot(5,2,8,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.02;
k=1;
normv = squeeze(vsigstimspont(4,[3 2 1 4],:));% ./ sum(vsigstimspont(2,[1 3:4],:),2));
randv = mean(squeeze(mean(vsigstimspont(4,5:end,:))));
for d = 1:4
	hold all;
	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),'o-','color',cm(d,:),'markersize',4);
end
plot([0 size(normv,1)+1], randv * [1 1],'k--','linewidth',2);
text(1.3,.22,'random','HorizontalAlignment','right','fontsize',8)
%title({'stim responses','of projections'},'fontweight','normal');
set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'stim-only','spont-only','behav-only','stim-behav'});
set(gca,'XTickLabelRotation',30);
xlim([.5 size(normv,1)+.5]);
ylabel({'signal-to-noise ratio'});
axis square;


i=i+1;
hs{i}=my_subplot(5,2,9,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.02;
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
hs{i}=my_subplot(5,2,10,[.9 .62]);
hs{i}.Position(2)=hs{i}.Position(2)-.015;
hs{i}.Position(1)=hs{i}.Position(1)-.05;
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
tstim=1:700;
r = corr(y(ttest),mg(ttest));
plot(mg(ttest(tstim)),'r','linewidth',1)
%plot(y(ttest(tstim))*as,'color',.6*[1 1 1],'linewidth',1);
plot(fg(:,ttest(tstim))'*a,'b','linewidth',0.5)
text(0,1.2,'gain','color','r','fontangle','normal','fontsize',8);
text(.25,1.2,sprintf('face prediction of gain\n          r=%1.2f',face_pred_gain_r2(1)),...
	'color','b','fontangle','normal','fontsize',8);
%text(.3,0.1,sprintf('stim-behav shared dim, r=%1.2f',r),'color',.4*[1 1 1],'fontangle','normal','fontsize',8);
box off;
axis tight;
axis off;
plot([0 60/1.2],[1 1]*min(mg(ttest(tstim))),'k','linewidth',2)
text(0,-.02,'1 min','fontangle','normal','fontsize',8);


%normv(2:3,:) = squeeze(fitmult(:,id(k),:)) .* normv(4,:);


tstr{1} = 'Variance of face PCs';
tstr{2} = {'Stim-behav shared','         subspace'};
tstr{3} = 'Top dimension';
tstr{5} = 'Direction decoding';
tstr{4} = 'Orientation decoding';
for k = 6:length(hs)
    tstr{k} = '';
end

% -------------- LETTERS
hp=.1;
hy=1.24;
deffont=8;
for j = [1:length(hs)]
   if j<length(hs)
        hp0=hp;
        hy0 = hy;
   else
	   hp0=0.04;
	   hy0=hy;
    end
    hpos = hs{j}.Position;
    lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
    axes('position', lpos);
    text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
    
    axes('position', [lpos(1)+.04 lpos(2)-.005 lpos(3) lpos(4)]);
    text(0,0, tstr{j},'fontangle','normal','fontsize',deffont);
    axis([0 1 0 1]);
    axis off;
end
%%
print(fullfile(matroot, 'suppOriStimSpont.pdf'),'-dpdf');

