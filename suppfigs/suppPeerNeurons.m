function suppPeerNeurons(matroot)


load(fullfile(matroot,'expv_behavior_neurons.mat'));
load(fullfile(matroot,'expv_timedelay_neurons.mat'));
load(fullfile(matroot,'PCApred.mat'));
try
ephys = load(fullfile(matroot,'ephys_peers.mat'));
catch
end

%%
close all;
default_figure([1 1 7.25 6.]);

%%

load('cdat.mat');
ndat = size(cdat,1);

% PC colors
nPCs = 4;
cpc = max(0,colormap('spring'));
cpc = cpc(1:32,:);
cpc = cpc(round(linspace(1,size(cpc,1),nPCs)),:);

xh = .55;
yh = xh;

clf;
% ------ PEER PRED DIAGRAM ------ %
i=0;
i=i+1;
hs{i}=my_subplot(3,4,1,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + 0;
axis off;
hp=hs{i}.Position;
hp(2)=hp(2)-.1;
hp(4)=hp(4)*1.3;
hp(1)=hp(1)-.01;
%hp(3) = hp(3) * 1.6;
%hp(2)=hp(2)+.04;
%hp(4)=hp(4)*.7;
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
text(.5,1.1,'neural PCs','HorizontalAlignment','center','fontsize',8);
text(1,.95,{'single neuron','activity'},'HorizontalAlignment','center','fontsize',8);
text(0,1.3,'peer prediction','fontweight','normal','fontangle','italic','fontsize',10);
% ---------- EPHYS --------%
i=i+1;
hs{i} = my_subplot(3,4,2,[xh*1.3 yh*1.3]);
hs{i}.Position(2) = hs{i}.Position(2) + .0;
hs{i}.Position(1) = hs{i}.Position(1) + .04;
try
nPC = 2.^[0:9];
[~,isort] = sort(max(ephys.expv_neurons,[],2),'descend');
for j = 1:size(ephys.expv_neurons,1)
    semilogx(nPC, ephys.expv_neurons(isort(j),:), 'color', ephys.colors(isort(j),:));
	hold all;
	text(1.05,1-(j-1)*.1,ephys.tgrps{isort(j)},'color',ephys.colors(isort(j),:),'fontangle','normal');
end
set(gca,'xtick',2.^[0:3:10]);
axis([1 512 0 .7]);
box off;
ylabel({'variance explained','(test set)'});
xlabel('# of PCs');
grid on;
grid minor;
grid minor;
axis square;
catch
end

% ----- PEER PRED NEURONS ---------%
i=i+1;
hs{i}=my_subplot(3,4,6,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) - .03;
ndims1 = 2.^[0:10];
for d = 1:ndat
	semilogx(ndims1,expv_neurons(:,d),'color',cdat(d,:),'linewidth',0.5);
	hold all;
end
semilogx(ndims1, mean(expv_neurons,2),'k','linewidth',2);
grid on;
grid minor;
grid minor;
set(gca,'xtick',2.^[0:4:12],'ytick',[0:.1:0.3]);
axis tight;
box off;
axis square;
ylim([0 .3]);
ylabel({'variance explained'});
xlabel({'# of PCs'});
title('peer prediction','fontweight','normal')
text(-2.2,1.2, {'Single neuron analyses','on two-photon data:'},'fontangle','italic','fontsize',10);


i=i+1;
hs{i} = my_subplot(3,4,7,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) - .03;
ev=[squeeze(expv_behavior(1,:,1:3)) squeeze(expv_behavior(3,:,7))'];
for j = 1:ndat
    plot([1:4],ev(j,:),'color',cdat(j,:))
    hold all;
end
plot([1:4],mean(ev),'ko-','linewidth',2,'markerfacecolor','k','markersize',4)
box off;
set(gca,'xtick',[1:8],'xticklabel',{'running','pupil area','whisking','all 3'});
xtickangle(45);
axis square;
axis tight;
ylabel({'variance explained'});
xlim([0.5 4.5]);
ylim([0 0.065]);


% --------------------- DIMS VS VAR EXP -------------- %
i=i+1;
hs{i} = my_subplot(3,4,8,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) - .03;
%hs{i}.Position(2) = hs{i-1}.Position(2);
for j = 1:ndat
    semilogx(ndims0,expv_behavior(:,j,9),'color',cdat(j,:));
    hold all;
end
box off;
set(gca,'xtick',[1 4 16 64]);
xlabel('dimensions');
ylabel({'variance explained'});
title('face prediction','fontweight','normal')
axis tight;
ylim([0 .15]);
xlim([0 64]);
hold all;
semilogx(ndims0,mean(expv_behavior(:,:,9),2),'k','linewidth',2);
axis square;
grid on;
grid minor;
grid minor;

%-------------- PEER PRED -----------------%
i=i+1;
hs{i}=my_subplot(3,4,9,[xh yh]);
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
hs{i} = my_subplot(3,4,10,[xh yh]);
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
hs{i} = my_subplot(3,4,11,[xh yh]);
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
hs{i} = my_subplot(3,4,12,[xh yh]);
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

%
% -------------- LETTERS
hp=.07;
hy=1.23;
deffont=8;
for j = [1:length(hs)]
    hp0=hp;
	hy0=hy;
	if j==2
		hy0 = 1.12;
	end
    hpos = hs{j}.Position;
    axes('position', [hpos(1)-hp0 hpos(2)+hpos(4)*hy0(1) .01 .01]);
    text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
end


%%
print(fullfile(matroot, 'suppPeerNeurons.pdf'),'-dpdf');
