load(fullfile(matroot,'example_behavior.mat'));
load(fullfile(matroot,'expv_behavior_PC.mat'));
load(fullfile(matroot,'PCpredict.mat'));

%%
default_figure([1 1 4 4]);

%%
tstr   = {'running','pupil area', 'whisking', 'running+pupil','running+whisker',...
    'pupil+whisker','running+pupil+whisker'};

% cm is running pupil whisker

cm(1,:) = [.2 .8 .2];
cm(2,:) = [.5 .6 .5];
cm(3,:) = [.1 .9 .6];
cpc = max(0,colormap('spring'));
cpc = cpc(1:32,:);
nPCs = 4;
cpc = cpc(round(linspace(1,size(cpc,1),nPCs)),:);
load('cdat.mat');
ndat = size(cdat,1);

% independent contributions
%indvar = squeeze(expv_behavior(3,:,7) - expv_behavior(2,:,4:6));
%indvar = indvar(:, [3 2 1]);

i = 0;

clf;
trange = 7300+[1:500];
i=i+1;
hs{i} = my_subplot(2,1,1,[.87 .55]);
%hs{i}.Position(3)=hp.Position(3);
xp = zscore(xex(trange,:),1,1);
xp = xp - min(xp,[],1);
xp(:,1) = xp(:,1)/1.3;
xp = min(5, xp);
%xp = xp ./ max(xp,[],1);
ht=text(-.07,0.07,'normalized units','fontsize',8,'fontweight','bold','fontangle','normal');
ht.Rotation=90;
ht=text(-.07,-1.55,'correlation','fontsize',8,'fontweight','bold','fontangle','normal');
ht.Rotation=90;
%ht=text(-.06,-5.1,'weights','fontsize',8,'fontweight','bold','fontangle','normal');
%ht.Rotation=90;
for j = [3 1 2]
    plot([0:1/3:1/3*(numel(trange)-1)],xp(:,j),'color',cm(j,:),'linewidth',1);
    hold all;
    text(.02,1.-(j-1)*.12,tstr{j},'color',cm(j,:),'fontweight','bold',...
        'fontsize',8,'fontangle','normal');
end
xlabel('time (s)');
axis tight;
box off;
set(gca,'color','none','YColor','none');
%ylabel('normalized units');

i=i+1;
hs{i} = my_subplot(2,2,3,[.6 .6]);
hs{i}.Position(2) = hs{i}.Position(2) + .02;
cmat = mean(cc,3);
cmat = cmat - diag(NaN*diag(cmat));
cmat = [squeeze(cc(1,2,:)) squeeze(cc(1,3,:)) squeeze(cc(2,3,:))];
for j = 1:ndat-1
    plot(cmat(j,:),'color',cdat(j,:))
    hold all;
end
plot(mean(cmat(1:end-1,:),1),'ko-','linewidth',2,'markerfacecolor','k','markersize',4)
axis tight;
box off;
ylim([0 .85]);
set(gca,'xtick',[1 2 3],'xticklabel',...
    {'run-pupil','run-whisk','pupil-whisk'});
ylabel('correlation');
axis square;
xtickangle(45);
xlim([0.75 3.25]);

% pcolor([cmat nan(3,1);nan(1,4)])
% shading flat;
% set(gca,'ydir','reverse');
% colormap(hs{i},'jet');
% set(gca,'xtick',[1 2 3]+.5,'xticklabel',{'running','pupil','whisking'});
% set(gca,'ytick',[1 2 3]+.5,'yticklabel',{'running','pupil','whisking'});
% colorbar;
% title('correlation');

evs = [];
id = [1 1 1 3 1];
b = [1 2 3 7 7];
for j = 1:5
	evs = [evs ((cov_neur(1,:)'-squeeze(cov_res_beh(1,id(j),:,b(j))))./var_neur(1,:)')];
end
	%std((cov_neur-squeeze(cov_res_beh(:,3,:,7)))./var_neur,1,2)/sqrt(ndat-1),{'color',[0 0.2 0],'linewidth',1})

i=i+1;
hs{i} = my_subplot(2,2,4,[.6 .6]);
hs{i}.Position(2) = hs{i}.Position(2) + .02;
%hs{i}.Position(3)=hs{i}.Position(3)+.03;
for j = 1:ndat
    plot([1:4],evs(j,1:4),'color',cdat(j,:))
    hold all;
end
plot([1:4],mean(evs(:,1:4)),'ko-','linewidth',2,'markerfacecolor','k','markersize',4)
box off;
set(gca,'xtick',[1:8],'xticklabel',{'running','pupil area','whisking','all 3'});
xtickangle(45);
axis square;
axis tight;
ylabel({'variance explained','(test set)'});
xlim([0.5 4.5]);
ylim([0 1]);
title('PC 1','fontsize',8,'fontweight','normal')


%
% -------------- LETTERS
hp=.08;
hy=1.23;
deffont=8;
for j = [1:length(hs)]
    if j==1
        hp0 = hp-.02;
    else
        hp0=hp;
    end
    hpos = hs{j}.Position;
    axes('position', [hpos(1)-hp0 hpos(2)+hpos(4)*hy(1) .01 .01]);
    text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
end

%%
print(fullfile(matroot,'supp1DBeh.pdf'),'-dpdf')