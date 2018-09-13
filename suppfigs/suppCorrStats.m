function suppCorrStats(matroot)

load(fullfile(matroot,'corr1stpc.mat'));

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
hs{i} = my_subplot(5,4,5,[xh yh]);
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
hs{i} = my_subplot(5,4,6,[xh yh]);
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
hs{i} = my_subplot(5,4,7,[xh yh]);
plot(results.cc1(1:50:end), results.cc2(1:50:end), '.','markersize',4,...
    'color',cdat(1,:));
axis tight;
box off;
xlabel('correlation half 1');
ylabel('correlation half 2');
text(.25,1.05, sprintf('r=%1.2f',results.rcorr(dex)));
axis square;

% ------------ R OF CORRELATIONS
i = i+1;
hs{i} = my_subplot(5,4,8,[xh yh]);
histogram(results.rcorr,'binedges',[.5:.1:.9],'edgecolor',[0 0 0],...
    'facecolor', .1*[1 1 1]);
box off;
ylabel('number of recordings');
xlabel('r of correlations');
xlim([0 1]);
axis square;

