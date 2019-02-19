function suppIncNeuronsBins(matroot)

%%

%%
close all;
default_figure([1 1 8 3.75]);
%%
clf;

xh=.65;
yh=.65;

i=0;


blu = [0 0 1];
red = [1 0 0];
green = [0 .5 0];

results=load(fullfile(matroot,'increasing_neurons_SVCA_facepred.mat'));
nneur0 = results.nneur0;

nl = length(nneur0);
col = [linspace(blu(1), red(1), nl); ...
	linspace(blu(2), red(2), nl); linspace(blu(3), red(3), nl)]';
cm = col;

i=i+1;
hs{i}=my_subplot(2,5,1, [xh yh]);
for j = 1:length(nneur0)
	semilogx(100*nanmean(nanmean(results.cov_neur(:,j,:,:)./results.var_neur(:,j,:,:),4),3),...
		'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
	text(1,.5+j*.1,sprintf('%d',nneur0(j)*2),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([0 100]);
xlim([1 1024]);
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]*100);
ylabel('% reliable variance');
xlabel('SVC dimension');
set(gca,'fontsize',6);
grid on;
grid minor;
grid minor;
axis square;
box off;

text(-.1,1.3,'Varying numbers of neurons (1.2 second time bins)','fontsize',10);
	
i=i+1;
hs{i}=my_subplot(2,5,2, [xh yh]);
for j = 1:length(nneur0)
	cv = nanmean(nanmean(results.cov_neur(:,j,:,:)./nansum(results.cov_neur(:,j,:,:),1),4),3);
	loglog(cv(1:min(nneur0(j),1024)),...
		'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
%	text(1,.5+j*.1,sprintf('%d',nneur0(j)*2),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([1e-5 1]);
xlim([1 1024]);
set(gca,'xtick',10.^[0:3],'ytick',10.^[-4:1:0]);
ylabel({'normalized','reliable variance'});
xlabel('SVC dimension');
grid on;
grid minor;
grid minor;
axis square;
box off;
set(gca,'fontsize',6);

i=i+1;
hs{i}=my_subplot(2,5,3, [xh yh]);
for j = 1:length(nneur0)
	cov_neur = squeeze(results.cov_neur(1:128,j,:,:));
	var_neur = squeeze(results.var_neur(1:128,j,:,:));
	plot(nneur0(j)*2,100*mean(mean(sum(cov_neur) ./ sum(var_neur))),'.', 'color', cm(j,:),'markersize',10);
	hold all;
	%text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
%ylim([0 100]);
%xlim([1 128]);
set(gca,'xtick',[2000:2000:1e4])%,'ytick',[0:1]*100);
set(gca,'fontsize',6);
ylabel('% variance explained');
xlabel('# of neurons');
title('SVC 1-128','fontweight','normal')
%ylim([0 20]);
axis tight;
ylim([0 70]);
grid on;
grid minor;
grid minor;
axis square;
box off;
%xlim([500 1e4]);

i=i+1;
hs{i}=my_subplot(2,5,4, [xh yh]);
for j = 1:length(nneur0)
	cov_neur = squeeze(results.cov_neur(1:128,j,:,:));
	var_neur = squeeze(results.var_neur(1:128,j,:,:));
	% put the lambda as the last index
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,5,j,:,:,:)), [1 3 4 2]);
	[~,ilam] = max(nanmean(nanmean(nansum((cov_neur - cov_res_beh)./var_neur,1),2),3));
	cov_res_beh = cov_res_beh(:,:,:,ilam);
	semilogx(100*nanmean(nanmean(results.cov_neur(:,j,:,:)./results.var_neur(:,j,:,:),4),3),'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
	semilogx(100*nanmean(nanmean((cov_neur - cov_res_beh)./var_neur, 2),3),'color', max(0,cm(j,:)-.2),'linewidth',1.5)
%	text(1,.5+j*.1,sprintf('%d',nneur0(j)*2),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([0 100]);
xlim([1 128]);
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]*100);
ylabel('% variance explained');
xlabel('SVC dimension');
title('from face motion','fontweight','normal')
grid on;
grid minor;
grid minor;
axis square;
box off;
set(gca,'fontsize',6);

i=i+1;
hs{i}=my_subplot(2,5,5, [xh yh]);
for j = 1:length(nneur0)
	cov_neur = squeeze(results.cov_neur(1:128,j,:,:));
	var_neur = squeeze(results.var_neur(1:128,j,:,:));
	% put the lambda as the last index
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,5,j,:,:,:)), [1 3 4 2]);
	[~,ilam] = max(nanmean(nanmean(nansum((cov_neur - cov_res_beh)./var_neur,1),2),3));
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,:,j,ilam,:,:)), [1 3 4 2]);
	semilogx(results.ndims0, 100*squeeze(nanmean(nanmean(nansum(cov_neur - cov_res_beh,1)./nansum(var_neur,1),2),3)),...
		'color', max(0,cm(j,:)-.2),'linewidth',1.5);
	hold all;
	%text(1,.5+j*.1,sprintf('%d',nneur0(j)*2),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
xlabel('rank');
ylabel({'% variance explained'});
title('SVC 1-128','FontWeight','normal');
axis tight;
ylim([0 23]);
xlim([0 128]);
set(gca,'xtick',[1 4 16 64]);
set(gca,'fontsize',6);
grid on;
grid minor;
grid minor;
axis square;
box off;


% ---------------- TIMEBINS CHANGE ----------- %


results=load(fullfile(matroot,'timebins_30Hz_spont.mat'));

tbins = results.tbins;
tbinss = tbins * 1/30;

nl = length(tbins);
col2 = [linspace(blu(1), green(1), nl); ...
	linspace(blu(2), green(2), nl); linspace(blu(3), green(3), nl)]';
%cm = col2(end:-1:1,:);
cm = col2;
%cm(end-1:end,

i=i+1;
hs{i}=my_subplot(2,5,6, [xh yh]);
for j = 1:length(tbins)
	semilogx(100*nanmean(results.cov_neur(:,j,:)./results.var_neur(:,j,:),3),...
		'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
	text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
text(.74,1,'timebin (s):','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([0 100]);
xlim([1 1024]);
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]*100);
ylabel('% reliable variance');
xlabel('SVC dimension');
set(gca,'fontsize',6);
grid on;
grid minor;
grid minor;
axis square;
box off;

text(-.1,1.35,'Varying time bin size (single plane imaging @ 30Hz, ~900 neurons)','fontsize',10);

i=i+1;
hs{i}=my_subplot(2,5,7, [xh yh]);
for j = 1:length(tbins)
	loglog(nanmean(results.cov_neur(:,j,:)./nansum(results.cov_neur(:,j,:),1),3),...
		'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
%	text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([1e-5 1]);
xlim([1 1024]);
set(gca,'xtick',10.^[0:3],'ytick',10.^[-4:1:0]);
ylabel({'normalized','reliable variance'});
xlabel('SVC dimension');
grid on;
grid minor;
grid minor;
axis square;
box off;
set(gca,'fontsize',6);

i=i+1;
tbs={};
hs{i}=my_subplot(2,5,8, [xh yh]);
for j = 1:length(tbins)
	cov_neur = squeeze(results.cov_neur(1:128,j,:));
	var_neur = squeeze(results.var_neur(1:128,j,:));
	plot(tbins(j)/30,100*mean(sum(cov_neur) ./ sum(var_neur)),'.', 'color', cm(j,:),'markersize',10);
	hold all;
	tbs{j} = sprintf('%2.2f',tbinss(j));
	%text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
%ylim([0 100]);
%xlim([1 128]);
axis tight;
set(gca,'xtick',[0:4]);
%set(gca,'xtick',tbins(1:2:end)/30,'xticklabel',tbs(1:2:end))%,'ytick',[0:1]*100);
set(gca,'fontsize',6);
ylabel('% variance explained');
xlabel('time bin size');
title('SVC 1-128','fontweight','normal')
ylim([0 40]);
xlim([0 4.3]);
grid on;
grid minor;
grid minor;
axis square;
box off;


i=i+1;
hs{i}=my_subplot(2,5,9, [xh yh]);
for j = 1:length(tbins)
	cov_neur = squeeze(results.cov_neur(1:128,j,:));
	var_neur = squeeze(results.var_neur(1:128,j,:));
	% put the lambda as the last index
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,end,j,:,:)), [1 3 2]);
	[~,ilam] = max(nanmean(nansum((cov_neur - cov_res_beh)./var_neur,1),2));
	cov_res_beh = cov_res_beh(:,:,ilam);
	semilogx(100*nanmean(cov_neur./var_neur,2),'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
	semilogx(100*nanmean((cov_neur - max(0,cov_res_beh))./var_neur, 2),'color', max(0,cm(j,:)),'linewidth',1.5)
	hold all;
	%text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([0 100]);
xlim([1 128]);
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]*100);
set(gca,'fontsize',6);
ylabel('% variance explained');
title('from face motion','fontweight','normal');
xlabel('SVC dimension');
grid on;
grid minor;
grid minor;
axis square;
box off;

i=i+1;
hs{i}=my_subplot(2,5,10, [xh yh]);
for j = 1:length(tbins)
	cov_neur = squeeze(results.cov_neur(1:128,j,:));
	var_neur = squeeze(results.var_neur(1:128,j,:));
	% put the lambda as the last index
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,5,j,:,:)), [1 3 2]);
	[~,ilam] = max(nanmean(nansum((cov_neur - cov_res_beh)./var_neur,1),2));
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,:,j,ilam,:)), [1 3 2]);
	semilogx(results.ndims0, 100*squeeze(nanmean(nansum(cov_neur - cov_res_beh,1)./nansum(var_neur,1),2)),...
		'color', max(0,cm(j,:)-.2),'linewidth',1.5);
	hold all;
	%text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
xlabel('rank');
ylabel({'% variance explained'});
title('SVC 1-128','FontWeight','normal');
set(gca,'fontsize',6);
axis tight;
ylim([0 23]);
xlim([0 64]);
set(gca,'xtick',[1 4 16 64]);
grid on;
grid minor;
grid minor;
axis square;
box off;



for k = 1:length(hs)
	axes('position',hs{k}.Position);
	axis off;
	%if k < 4
		text(-.32,1.08,char(64+k),'fontsize',10,'fontweight','bold','fontangle','normal');
	%else
	%	text(-.4,1.05,char(64+k),'fontsize',10,'fontweight','bold','fontangle','normal');
	%end
	%text(-.2,1.15,tstr{j},'fontsize',8);
	
end


%%
print(fullfile(matroot, 'suppIncNeuronsBins.pdf'),'-dpdf');


