function suppCCEphys(matroot)

load(fullfile(matroot, 'ephysResults_KS2.mat'));

%% cross correlations
close all;
default_figure([0 0 4.25 4.25]);

%%
tlag = results.tlag;
id = [1:length(tlag)];
pneur=[];
pbeh=[];
for d=1:3
	cov_neur0 = squeeze(results.cov_neur_t(1:128,ceil(length(tlag)/2),d));
	var_neur0 = squeeze(results.var_neur(1:128,d));
	for j = 1:length(tlag)	
		
		cov_res_beh0 = squeeze(results.cov_res_beh_t(:,6,j,:,d));
		[m(j),ilam] = max(nansum(cov_neur0 - cov_res_beh0)./nansum(var_neur0));
		pbeh(j,d) = nansum(cov_neur0 - cov_res_beh0(:,ilam)) / nansum(var_neur0);
	end
	pneur(:,d) = sum(results.cov_neur_t(1:128, id, d)) ./ sum(var_neur0, 1);
end


pbeh = pbeh*100;
pneur=pneur*100;
%
xh=.6;
yh=xh;
clf
i=0;
i=i+1;
hs{i}=my_subplot(2,3,1,[xh yh]);
plot(tlag(id),pneur,'r');
xlabel('time lag (s)');
ylabel('% variance explained');
title('prediction from neurons','fontweight','normal','color','r');
axis tight;
box off;
axis square;
grid on;
set(gca,'fontsize',6);
ylim([0 40]);

i=i+1;
hs{i}=my_subplot(2,3,2,[xh yh]);
plot(tlag,pbeh,'b');
title('from face motion','fontweight','normal','color','b');
axis tight;
grid on;
ylim([0 8]);
box off;
axis square;
xlabel('time from behav. (s)');
ylabel('% variance explained');
%xlim([-2 2]);
grid on;
set(gca,'fontsize',6);

% i=i+1;
% hs{i}=my_subplot(1,5,3,[xh yh]);
% plot(tlag,mean(pbeh,2),'color',[0 0 .8]);
% hold all;
% plot(tlag,mean(pneur,2),'color',[.8 0 0]);
% box off;
% axis square;
% axis tight;
% ylim([0 10]);
% xlabel('time lag (s)');
% ylabel('% variance explained');
% title('averaged','fontweight','normal');

i=i+1;
hs{i}=my_subplot(2,3,3,[xh yh]);
plot(tlag,mean(pbeh,2) / max(mean(pbeh,2)),'color',[0 0 .8]);
hold all;

plot(tlag(id),mean(pneur,2) / max(mean(pneur,2)),'color',[.8 0 0]);
box off;
axis square;
axis tight;
ylim([0 1]);
xlabel('time lag (s)');
ylabel({'normalized','variance explained'});
title('normalized','fontweight','normal');
grid on;
set(gca,'fontsize',6);

i=i+1;
hs{i}=my_subplot(2,1,2);
plot(tlag,mean(pbeh,2) / max(mean(pbeh,2)),'color',[0 0 .8]);
hold all;
plot(tlag(id),mean(pneur,2) / max(mean(pneur,2)),'color',[.8 0 0]);
box off;
set(gca,'fontsize',6);
axis tight;
ylim([0 1]);
xlabel('time lag (s)');
ylabel({'normalized','variance explained'});
title('zoom in','fontweight','normal');
xlim([-0.8 0.8]);
grid on;
	hp=.08;
hy=1.1;
deffont=8;



for j=1:length(hs)
	hp0=hp; hy0=hy;
    if j==4
        hp0 = 0.04;
        hy0 = 1.15;
	end
    hpos = hs{j}.Position;
    lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
    axes('position', lpos);
    t=text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','FontName', 'Helvetica');
    axis([0 1 0 1]);
    axis off;
end

%%
print(fullfile(matroot,'suppTimelags.pdf'),'-dpdf')

%%

for k = 1:2
	for d = 1:3
if k==1
	tl = pneur(:,d);
else
	tl = pbeh(:,d);
end

tli = interp1(tlag(:),tl(:),[-8:.01:8]');

[~,ix1]=min(abs(max(tli,[],1)/2 - tli(1:800,:)));
[~,ix2]=min(abs(max(tli,[],1)/2 - tli(801:end,:)));

fwhm(k,d) = (800-ix1+ix2)*.01;

	end
end
fwhm
