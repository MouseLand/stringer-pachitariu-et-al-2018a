function suppCorrStats(matroot)

%%
load(fullfile(matroot,'corr1stpc.mat'));

%%
close all;
default_figure([1 1 6.25 1.75]);
%%
clf;
load('cdat.mat');
xh = 0.5; yh = 0.5;
i=0;
%
% -------------------- CORR MATRICES ---------%
ineu1 = 1:20:20*10;
ineu2 = 2:20:20*10;
rng(6);
ineu1 = randperm(size(results.cc1,1),10);
ineu1 = randperm(size(results.cc1,1),10);
crange = [-.15 .15];
i = i+1;
hs{i} = my_subplot(1,4,1,[xh yh]);
axis off;
hp=hs{i}.Position;
hp(1)=hp(1)-.03;
hp(2)=hp(2)-.01;
hp(3)=hp(3)+.025;
hp(4)=hp(4)+.02;
axes('position',hp);
imagesc(results.cc1(ineu1,ineu2), crange);
text(.5, 0,'neurons 11-20','HorizontalAlignment','center','fontangle','normal','fontsize',8);
ht=text(-.17, .5,'neurons 1-10','HorizontalAlignment','center','fontangle','normal','fontsize',8);
ht.Rotation=90;
axis off;
axis square;
colormap(gca,'redblue');
axes('position',[hp(1)+hp(3) hp(2) .03 hp(4)]);
cb=colorbar;
cb.Ticks=[1/6 .5 5/6];
cb.TickLabels={'-0.1','0','0.1'};
colormap(gca,'redblue');
axis off;

i = i+1;
hs{i} = my_subplot(1,4,2,[xh yh]);
axis off;
hp=hs{i}.Position;
hp(1)=hp(1)-.03;
hp(2)=hp(2)-.01;
hp(3)=hp(3)+.025;
hp(4)=hp(4)+.02;
axes('position',hp);
imagesc(results.cc2(ineu1,ineu2), crange);
axis off;
axis square;
colormap(gca,'redblue');
text(.5, 0,'neurons 11-20','HorizontalAlignment','center','fontangle','normal','fontsize',8);
ht=text(-.17, .5,'neurons 1-10','HorizontalAlignment','center','fontangle','normal','fontsize',8);
ht.Rotation=90;

% --- corr halves X Y ---- %
i = i+1;
hs{i} = my_subplot(1,4,3,[xh yh]);
plot(results.cc1(1:50:end), results.cc2(1:50:end), '.','markersize',4,...
    'color',cdat(dex,:));
axis tight;
box off;
xlabel('correlation half 1');
ylabel('correlation half 2');
text(.25,1.05, sprintf('r=%1.2f',results.rcorr(dex)),'fontangle','normal');
axis square;

% ------------ R OF CORRELATIONS
i = i+1;
hs{i} = my_subplot(1,4,4,[xh yh]);
histogram(results.rcorr,'binedges',[.5:.1:.9],'edgecolor',[0 0 0],...
    'facecolor', .1*[1 1 1]);
box off;
ylabel('# of recordings');
xlabel('r of correlations');
xlim([0 1]);
axis square;

%%
tstr{1} = 'correlation half 1';
tstr{2} = 'correlation half 2';
tstr{3} = '';
tstr{4} = '';

% -------------- LETTERS
hp=.06;
hy=1.25;
deffont=8;
for j = [1:length(hs)]
    hp0=hp;
    hy0 = hy;
    
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
print(fullfile(matroot,'suppCorrStats.pdf'),'-dpdf')