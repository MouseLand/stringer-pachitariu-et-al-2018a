function fig4new(matroot)

load(fullfile(matroot,'stimvar_30Hz.mat'));
load(fullfile(matroot,'faceSpectrum_30Hz.mat'));

%%
close all;
default_figure([8 1 8 5]);

%%
rng('default');
ktype = 1; % behavior space

clear hs;
i=0;
clf;
cm = colormap('hot');
cm = cm(32+[1:4:4*4],:);

dx = .55;
dy = .55;

i=i+1;
hs{i} = my_subplot(3,4,1,[dx dy]);
hold all;
for d = 1:size(Vshared,2)
    plot(Vshared(1:10,d,ktype,2)*100,'.-','markersize',8,'color',cm(d,:),'linewidth',.25)
	plot(Vshared(1:10,d,ktype,1)*100,'x-','markersize',8,'color',cm(d,:),'linewidth',.25)
end
text(1,1,'x = 200 ms','fontsize',8);
box off;
ylabel({'% stim variance'},'fontsize',8);
xlabel('dimension','fontsize',8);
%title('stim-face subspace','fontweight','normal','fontangle','italic','fontsize',10)
axis square;
axis tight;
ylim([0 8]);
set(gca,'xtick',[1 5 10]);


i=i+1;
hs{i} = my_subplot(3,4,2,[dx dy]);
hold all;
for d = 1:size(Ushared,1)
	for tb = 1:2
    uplot = Ushared{d,ktype,tb};
    uplot = uplot * sign(mean(uplot));
    bed = [-.02:.001:.035];
    histogram(uplot,20,'binedges',bed,'edgecolor',cm(d,:),...
        'normalization','probability','displaystyle','stairs');
	
	end
end
box off;
axis tight;
xlim([-.01 .03]);
ylim([0 0.2]);
ylabel('fraction','fontsize',8);
xlabel('neural weights','fontsize',8);
%title('1st PC weights','fontweight','normal','fontangle','italic','fontsize',10)
axis square;

i=i+1;
hs{i}=my_subplot(3,4,3,[dx dy]);
k=1;
tm = {'x-','o-'};
for tb=1:2
normv = squeeze(decoding([3 2 1],:,1,tb));% ./ sum(vsigstimspont(2,[1 3:4],:),2));
randv = mean(decoding(5:end,:,1),1,tb);
normv = [normv(1:3, :); randv];% normv(end,:)];
normv = normv * 360/33;
for d = 1:4
	hold all;
	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),tm{tb},'color',cm(d,:),'markersize',4);
end
end
%title({'stim responses','of projections'},'fontweight','normal');
set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'stim-only','spont-only','behav-only','random 32D'});
set(gca,'XTickLabelRotation',30);
xlim([.5 size(normv,1)+.5]);
ylabel({'decoding error'});
axis square;

i=i+1;
hs{i}=my_subplot(3,4,5,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.025;
cp=[.4 .6 1; 0 0 0; 0.8 0.3 0.8; .5 .5 .5];
tm = {'x','o'};
for tb=1:2
ll=[squeeze(vsigstimspont(3,:,:,tb)),squeeze(vsigstimspont(2,:,:,tb))];
ll = ll(:);
loglog([min(ll) max(ll)],[min(ll) max(ll)],'k','linewidth',.5)
hold all;
for k=1:4
	loglog(squeeze(vsigstimspont(3,k,:)),squeeze(vsigstimspont(2,k,:)),tm{tb},'color',cp(k,:))
	hold all;
end
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
hs{i}=my_subplot(3,4,6,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.025;
hold all;
k=0;
tb = 1;
for j = isti
    k=k+1;
    plot(projstim{1}{3}{tb}(istims{1}==j,isti(1)),projstim{1}{3}{tb}(istims{1}==j,isti(2)),...
		'.','color',cs(k,:),'markersize',5);
	plot(Rfit{1}{3}{tb}(istims{1}==j,isti(1)), Rfit{1}{3}{tb}(istims{1}==j,isti(2)), 'r','linewidth',0.5);
    %text(1.1,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
end
text(1.2,1,'stim 1', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.2,0.85,'stim 2', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.2,0.7,{'multiplicative','gain model'},'color','r','HorizontalAlignment','right','fontangle','normal','fontsize',8);
%axis tight;
axis([-80 180 -80 200]);
%axis([0 220 0 250]);
axis square;
box off;
title('200 ms bins');
xlabel('stim-only proj1');
ylabel('stim-only proj2');

i=i+1;
hs{i}=my_subplot(3,4,7,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.02;
hold all;
k=0;
for j = isti
    k=k+1;
    plot(projstim{1}{ktype}{tb}(istims{1}==j,1),projstim{1}{ktype}{tb}(istims{1}==j,2)*-1,...
		'.','color',cs(k,:),'markersize',5);
	%plot(Rfit{1}{1}(istims{1}==j,isti(1)), Rfit{1}{1}(istims{1}==j,isti(2)), 'k','linewidth',0.5);
    %text(1.3,1-(k-1)*.15,sprintf('stim %d',k),'color',cs(k,:),'fontangle','normal','fontsize',8,'HorizontalAlignment','right')
end
axis tight;
axis square;
box off;
title('200 ms bins');
xlabel('behav-only proj1');
ylabel('behav-only proj2');
text(1.4,1,'stim 1', 'color',cs(1,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);
text(1.4,0.85,'stim 2', 'color',cs(2,:),'HorizontalAlignment','right','fontangle','normal','fontsize',8);


i=i+1;
hs{i}=my_subplot(3,4,8,[dx dy]);
%hs{i}.Position(2)=hs{i}.Position(2)+0.02;
k=1;
tm = {'x-','o-'};
for tb = 1:2
normv = squeeze(vsigstimspont(4,[3 2 1 4],:,tb));% ./ sum(vsigstimspont(2,[1 3:4],:),2));
randv = mean(squeeze(mean(vsigstimspont(4,5:end,:,tb))));
for d = 1:4
	hold all;
	plot([1:size(normv,1)]+.05*randn(1,size(normv,1)), normv(:,d),tm{tb},'color',cm(d,:),'markersize',4);
end
plot([0 size(normv,1)+1], randv * [1 1],'k--','linewidth',2);
end
text(1.3,.22,'random','HorizontalAlignment','right','fontsize',8)
%title({'stim responses','of projections'},'fontweight','normal');
set(gca,'xtick',[1:size(normv,1)],'xticklabel',{'stim-only','spont-only','behav-only','stim-behav'});
set(gca,'XTickLabelRotation',30);
xlim([.5 size(normv,1)+.5]);
ylabel({'signal-to-noise ratio'});
axis square;


i=i+1;
hs{i}=my_subplot(3,4,9,[dx dy]);
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
hs{i}=my_subplot(3,4,10,[.9 .62]);
%hs{i}.Position(2)=hs{i}.Position(2)-.015;
hs{i}.Position(1)=hs{i-4}.Position(1);
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
tstr{4} = '';
for k = 5:length(hs)
    tstr{k} = '';
end

% -------------- LETTERS
hp=.05;
hy=1.24;
deffont=8;
for j = [1:length(hs)]
   if j<length(hs)
        hp0=hp;
        hy0 = hy;
   else
	   hp0=0.03;
	   hy0=hy;
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
print(fullfile(matroot, 'suppOriStimSpont.pdf'),'-dpdf');

