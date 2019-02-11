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

cm = colormap('hsv');
cm = cm(1:10:end,:);
results=load(fullfile(matroot,'increasing_neurons_SVCA_facepred.mat'));
nneur0 = results.nneur0;

i=i+1;
hs{i}=my_subplot(2,4,1, [xh yh]);
for j = 1:length(nneur0)
	semilogx(nanmean(nanmean(results.cov_neur(:,j,:,:)./results.var_neur(:,j,:,:),4),3),...
		'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
	text(1,.5+j*.1,sprintf('%d',nneur0(j)*2),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([0 1]);
xlim([1 1024]);
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]);
ylabel('% reliable variance');
xlabel('SVC dimension');
grid on;
grid minor;
grid minor;
axis square;
box off;

text(-.1,1.2,'Varying numbers of neurons (1.2 second time bins)','fontsize',10);
	
i=i+1;
hs{i}=my_subplot(2,4,2, [xh yh]);
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

i=i+1;
hs{i}=my_subplot(2,4,3, [xh yh]);
for j = 1:length(nneur0)
	cov_neur = squeeze(results.cov_neur(1:128,j,:,:));
	var_neur = squeeze(results.var_neur(1:128,j,:,:));
	% put the lambda as the last index
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,5,j,:,:,:)), [1 3 4 2]);
	[~,ilam] = max(nanmean(nanmean(nansum((cov_neur - cov_res_beh)./var_neur,1),2),3));
	cov_res_beh = cov_res_beh(:,:,:,ilam);
	semilogx(nanmean(nanmean(results.cov_neur(:,j,:,:)./results.var_neur(:,j,:,:),4),3),'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
	semilogx(nanmean(nanmean((cov_neur - cov_res_beh)./var_neur, 2),3),'color', max(0,cm(j,:)-.2),'linewidth',1.5)
%	text(1,.5+j*.1,sprintf('%d',nneur0(j)*2),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([0 1]);
xlim([1 128]);
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]);
ylabel('variance explained');
xlabel('SVC dimension');
grid on;
grid minor;
grid minor;
axis square;
box off;

i=i+1;
hs{i}=my_subplot(2,4,4, [xh yh]);
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
grid on;
grid minor;
grid minor;
axis square;
box off;


% ---------------- TIMEBINS CHANGE ----------- %

cm = colormap('hsv');
cm = cm(1:8:end,:);
results=load(fullfile(matroot,'timebins_30Hz_spont.mat'));
ibin = [1 2 3 5 6 8];%1:8;

tbins = results.tbins(ibin);
tbinss = tbins * 1/30;
results.cov_neur = results.cov_neur(:,ibin,:);
results.var_neur = results.var_neur(:,ibin,:);
results.cov_res_beh = results.cov_res_beh(:,:,ibin,:,:);

i=i+1;
hs{i}=my_subplot(2,4,5, [xh yh]);
for j = 1:length(tbins)
	semilogx(nanmean(results.cov_neur(:,j,:)./results.var_neur(:,j,:),3),...
		'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
	text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
text(.74,1,'timebin (s):','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([0 1]);
xlim([1 1024]);
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]);
ylabel('% reliable variance');
xlabel('SVC dimension');
grid on;
grid minor;
grid minor;
axis square;
box off;

text(-.1,1.2,'Varying time bin size (single plane imaging @ 30Hz, ~900 neurons)','fontsize',10);

i=i+1;
hs{i}=my_subplot(2,4,6, [xh yh]);
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

i=i+1;
hs{i}=my_subplot(2,4,7, [xh yh]);
for j = 1:length(tbins)
	cov_neur = squeeze(results.cov_neur(1:128,j,:));
	var_neur = squeeze(results.var_neur(1:128,j,:));
	% put the lambda as the last index
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,5,j,:,:)), [1 3 2]);
	[~,ilam] = max(nanmean(nansum((cov_neur - cov_res_beh)./var_neur,1),2));
	cov_res_beh = cov_res_beh(:,:,ilam);
	semilogx(nanmean(cov_neur./var_neur,2),'color', min(1,cm(j,:)+.4),'linewidth',0.5)
	hold all;
	semilogx(nanmean((cov_neur - max(0,cov_res_beh))./var_neur, 2),'color', max(0,cm(j,:)),'linewidth',1.5)
	hold all;
	%text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
ylim([0 1]);
xlim([1 128]);
set(gca,'xtick',10.^[0:3],'ytick',[0:.2:1]);
ylabel('variance explained');
xlabel('SVC dimension');
grid on;
grid minor;
grid minor;
axis square;
box off;

i=i+1;
hs{i}=my_subplot(2,4,8, [xh yh]);
for j = 1:length(tbins)
	cov_neur = squeeze(results.cov_neur(1:128,j,:));
	var_neur = squeeze(results.var_neur(1:128,j,:));
	% put the lambda as the last index
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,5,j,:,:)), [1 3 2]);
	[~,ilam] = max(nanmean(nansum((cov_neur - cov_res_beh)./var_neur,1),2));
	cov_res_beh = permute(squeeze(results.cov_res_beh(:,:,j,ilam,:)), [1 3 2]);
	loglog(results.ndims0, 100*squeeze(nanmean(nansum(cov_neur - cov_res_beh,1)./nansum(var_neur,1),2)),...
		'color', max(0,cm(j,:)-.2),'linewidth',1.5);
	hold all;
	%text(1,.4+j*.1,sprintf('%2.2f',tbinss(j)),'color',cm(j,:),'fontsize',8,'fontangle','normal','HorizontalAlignment','right')
end
%text(.74,1,'# of neurons:','fontsize',8,'fontangle','normal','HorizontalAlignment','right')
xlabel('rank');
ylabel({'% variance explained'});
title('SVC 1-128','FontWeight','normal');
axis tight;
ylim([0 23]);
xlim([0 64]);
set(gca,'xtick',[1 4 16 64],'ytick',2.^[-1:6]);
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


