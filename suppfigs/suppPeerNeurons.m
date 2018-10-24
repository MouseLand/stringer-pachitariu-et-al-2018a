
% ------ PEER PRED DIAGRAM ------ %
i=i+1;
hs{i}=my_subplot(5,5,23,[xh yh]);
axis off;
hold all;
msize = 6;
nls = [6 3 6];
call{1} = cpc;
call{2} = repmat([0 0 0],nls(2),1);
call{3} = repmat([1 1 1]*.8,nls(3),1);
ells = [3 3 3];
hp=hs{i}.Position;
hp(1)=hp(1)-.02;
hp(3)=hp(3)+.03;
hp(2)=hp(2)-.04;
hp(4)=hp(4)+.02;
layeredNet(hp, nls, call, msize,ells);
text(0,1.23,'PCs','HorizontalAlignment','center','fontsize',8);
text(.5,1.03,'rank','HorizontalAlignment','center','fontsize',8);
text(1,1.23,{'prediction','',''},'HorizontalAlignment','center','fontsize',8);
axis off;

% ----- PEER PRED NEURONS ---------%
i=i+1;
hs{i}=my_subplot(5,5,24,[xh yh]);
hs{i}.Position(1)=hs{i}.Position(1)+.03;
for d = 1:ndat
	semilogx(ndims1,expv128(:,d),'color',cdat(d,:),'linewidth',0.5);
	hold all;
end
semilogx(ndims1, mean(expv128,2),'k','linewidth',2);
grid on;
grid minor;
grid minor;
set(gca,'xtick',2.^[0:2:12],'ytick',[0:.2:1]);
axis tight;
box off;
axis square;
ylim([0 1]);
ylabel({'explainable variance'});
xlabel({'regression rank'});
title('PC 1-128');

%%
%%


% ----- PEER PRED TOTAL ---------%
i=i+1;
hs{i}=my_subplot(5,4,20,[xh yh]);
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




% ------ PEER PRED DIAGRAM ------ %
i=i+1;
hs{i}=my_subplot(5,4,17,[xh yh]);
axis off;
hold all;
msize = 6;
nls = [8 4 1];
for j = 1:length(nls)
    call{j} = repmat([0 0 0],nls(j),1);
end
call{2} = cpc;
ells = [3 3 0];
hp=hs{i}.Position;
hp(1)=hp(1)-.05;
hp(3)=hp(3)+.03;
hp(2)=hp(2)-.04;
hp(4)=hp(4)+.02;
layeredNet(hp, nls, call, msize,ells);
text(0,1.23,'peers','HorizontalAlignment','center','fontsize',8);
text(.5,1.03,'neural PCs','HorizontalAlignment','center','fontsize',8);
text(1,0.5,{'single','neuron','activity'},'HorizontalAlignment','center','fontsize',8);
axis off;
