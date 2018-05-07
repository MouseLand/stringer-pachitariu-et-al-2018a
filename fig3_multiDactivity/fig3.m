function fig3(matroot)

try
    dd = load(fullfile(matroot,'PCApred.mat'));
catch
    dd = load('PCApred.mat');
end
load(fullfile(matroot,'clust1D.mat'));
expv_neurons = dd.expv_neurons;
ndims0 = 2.^[0:size(expv_neurons,1)-1];
beh = results.beh;

%%

close all;
default_figure([1 1 9 5]);


%%
i= 0;
clf;
clear hs;
% group colors
nC = size(results.depth,1);
cg = colormap('jet');
rng(6);
cg = cg(randperm(64),:);
cg = cg(round(linspace(1,64,nC)),:);

% dataset colors
load('cdat.mat');
ndat = size(cdat,1);

% PC colors
nPCs = 4;
cpc = max(0,colormap('spring'));
cpc = cpc(1:32,:);
cpc = cpc(round(linspace(1,size(cpc,1),nPCs)),:);

i=i+1;
x0 = .03;
y0 = .56;
xh = .96;
yh = .44;
hs{i} = axes('position',[x0 y0 xh yh]);
trange = 2900+[1:3000];
imagesc(zscore(my_conv2(results.spks(:,trange), [5 0.5], [1 2]),1,2),[0.1 0.7]);
hold all;
NN=size(results.spks,1);
plot([1 1]*-15, NN-1000+[0 1000],'k');
ht=text(-.023,-.03,'  1000 neurons','fontangle','normal','fontsize',8);
ht.Rotation = 90;
axis tight;
cmap=colormap('gray');
colormap(flipud(cmap));
axis off;
axes('position',[x0 y0-.21 xh .2]);
hold all;
cm(1,:) = [.2 .8 .2];
cm(2,:) = [.5 .6 .5];
cm(3,:) = [0 .2 0];
tstr = {'running','pupil area','whisking'};
beh(:,1)=min(3.5,beh(:,1));
sbeh = sign(skewness(beh));
p=[0 .5 0];
for j = 1:3
    plot(beh(trange,j)*sbeh(j)-2.75*(j-1)+p(j),'linewidth',1,'color',cm(j,:));
    text(.07+(j-1)*.1,1,tstr{j},'fontangle','normal','fontweight','bold','fontsize',8,'color',cm(j,:));
end
plot([0 60/1.2],[-1 -1]*3*j+2,'k')
text(0,0,'  1 min','fontangle','normal','fontsize',8);

axis off;
axis tight;


yd = -.12;
xh=.63;
yh=.63;

i=i+1;
hs{i}=my_subplot(2,6,7,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + yd;
axis off;
hp=hs{i}.Position;
axes('Position',[hp(1)-.01 hp(2) hp(3:4)*1.12]);
hold all;
iplane = 5;
for j = 1:nC
    plot(results.pos(results.iclust==j & results.pos(:,3)==iplane*35,1),...
        results.pos(results.iclust==j & results.pos(:,3)==iplane*35,2),'.',...
        'color',cg(j,:),'markersize',4);
end
plot([0 500],[0 0]-30,'k','linewidth',2)
text(0,0,'500 \mum','fontsize',8,'fontangle','normal');
axis tight;
box off;
axis square;
axis off;
%text(-.24,1.35, 'example plane','fontangle','normal');

i=i+1;
hs{i}=my_subplot(2,6,8,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + yd;
hold all;
plot([0 530],[0 530],'k');
for j = 1:nC
    plot([results.indist(j,dex) results.outdist(j,dex)],'color',cg(j,:));
end
axis([.75 2.25 0 550])
box off;
axis square;
ylabel('pairwise distance (\mum)')
set(gca,'xtick',[1 2],'xticklabel',{'in-group','out-of-group'});

i=i+1;
hs{i}=my_subplot(2,6,9,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + yd;
hold all;
for d = 1:ndat
    plot(mean([results.indist(:,d) results.outdist(:,d)],1),'color',cdat(d,:),'linewidth',.5);
end
axis([.75 2.25 0 550])
box off;
axis square;
ylabel('pairwise distance (\mum)')
set(gca,'xtick',[1 2],'xticklabel',{'in-group','out-of-group'});

i=i+1;
hs{i}=my_subplot(2,6,10,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + yd;
axis off;
hp=hs{i}.Position;
hp(1)=hp(1)-.01;
hp(3) = hp(3) * 1.6;
hp(2)=hp(2)+.04;
hp(4)=hp(4)*.7;
axes('position',hp);
axis off;
msize = 6;
nls = [8 4 1];
for j = 1:length(nls)
    call{j} = repmat([0 0 0],nls(j),1);
end
call{2} = cpc;
ells = [3 3 0];
layeredNet(hp, nls, call, msize,ells);
text(0,1.16,'peers','HorizontalAlignment','center','fontsize',8);
text(.5,1.03,'neural PCs','HorizontalAlignment','center','fontsize',8);
text(1,1,{'single neuron','activity'},'HorizontalAlignment','center','fontsize',8);

i=i+1;
hs{i}=my_subplot(2,6,12,[xh yh]);
axis off;
%hs{i}.Position(1) = hs{i}.Position(1) - .01;
hs{i}.Position(1) = hs{i}.Position(1) - .05;
hs{i}.Position(2) = hs{i}.Position(2) + yd;
hp=hs{i}.Position;
axes('position',[hp(1:2) hp(3:4)*1.33]);
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


tstr{1}='';
tstr{2}='';
tstr{3}='Per group';
tstr{4}='Averaged';
tstr{5}='Peer prediction';
tstr{6}={''};

% -------------- LETTERS
hp=.05;
hy=1.05;
deffont=8;
for j = [1:length(hs)]
    if j==1
        hp0=.03;
        hy0=1.02;
    elseif j==2
        hp0=0.04;
        hy0=hy;
    elseif j==length(hs)
        hp0=.07;
        hy0=hy;
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
    
    axes('position', [lpos(1)+.03 lpos(2)-.008 lpos(3) lpos(4)]);
    text(0,0, tstr{j},'fontangle','normal','fontsize',8);
    axis([0 1 0 1]);
    axis off;
    
end
