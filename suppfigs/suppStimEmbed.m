function suppStimEmbed(matroot)

load(fullfile(matroot,'faceSpectrum.mat'));
load(fullfile(matroot,'stimvar.mat'));
load(fullfile(matroot,'stimfaceRepresentation.mat'));
load(fullfile(matroot,'exResponses.mat'));

%%
close all;
default_figure([1 1 8 8/3]);
%%
clf;
dex = 1;
[sval,sbins] = sort(sv{dex});
[fval,fbins] = sort(facevar{dex});

xh= 0.55;
yh= 0.55;
i = 0;
i = i+1;
hs{i}=my_subplot(1,3,1,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + .04;
plot(sbins(1:3:end),fbins(1:3:end),'.','markersize',2);
xlabel({'from stimulus (rank)'});
ylabel({'from behavior (rank)'});
box off;
axis tight;
axis([1 max(sbins) 1 max(sbins)]);
axis square;
set(gca,'ytick',10.^[0 4],'yticklabel',{'1','10^4'},...
	'xtick',10.^[0 4],'xticklabel',{'1','10^4'});

i = i+1;
hs{i}=my_subplot(1,3,2,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + .04;
axis square;
axis off;
hp=hs{i}.Position;
hp(3)=.05;
hp(1)=hp(1)-.06;
hss=axes('position',hp);
ss = sface{dex}(iface{dex},:);
ss = my_conv2(ss,10,1);
ss = my_conv2(ss,2,1);
imagesc(ss,[-.5 .5])
set(gca,'ytick',10.^[0 4],'yticklabel',{'1','10^4'},...
	'xtick',[]);
box off;
text(.5,-.05,{'behavior','coefficients'},'HorizontalAlignment','center','fontsize',8,'fontangle','normal');
ylabel('behavior embedding rank');
colormap(hss,redblue);

hp2=hp;
hp2(1)=hp2(1)+hp(3)+.04;
hp2(3)=hs{i}.Position(3);
%hp2(1)=hp(1);
%hp2(3)=hp2(4);
axes('position',hp2);
plot(istim{dex}(1:3:end),iface{dex}(1:3:end),'.','markersize',2);
box off;
axis tight;
axis([1 max(sbins) 1 max(sbins)]);
axis square;
set(gca,'ytick',10.^[0 4],'yticklabel',{'1','10^4'},...
	'xtick',10.^[0 4],'xticklabel',{'1','10^4'},...
	'ydir','reverse');

hp3=hp2;
hp3(4)=.13;
hp3(2)=hp3(2)-hp3(4)-.075;
axes('position',hp3);
ss = sstim{dex}(istim{dex},:);
ss = my_conv2(ss,10,1);
[~,imax] = max(ss,[],1);
[~,isort] = sort(imax);
ss = sstim{dex}(istim{dex},:);
ss = my_conv2(ss,3,1);
imagesc(ss(:,isort(end:-1:1))',[-.5 .5])
set(gca,'xtick',10.^[0 4],'xticklabel',{'1','10^4'},...
	'ytick',[]);
box off;
ht=text(-.11,.5,'stimuli','HorizontalAlignment','center','fontsize',8,'fontangle','normal');
ht.Rotation=90;
xlabel('stimulus embedding rank');
colormap(gca,redblue);


i = i+1;
hs{i}=my_subplot(1,3,3,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) + .04;
pds = abs(istim{dex}(1:7:end) - istim{dex}(1:7:end)');
pds = pds - diag(NaN*diag(pds));
pdf = abs(iface{dex}(1:7:end) - iface{dex}(1:7:end)');
pdf = pdf - diag(NaN*diag(pdf));
plot(pds(1:501:end), pdf(1:501:end), '.','markersize',2);
set(gca,'ytick',10.^[0 4],'yticklabel',{'1','10^4'},...
	'xtick',10.^[0 4],'xticklabel',{'1','10^4'});
p = polyfit(pds(~isnan(pds)),pdf(~isnan(pds)),1);
hold all;
plot([0 1.2e4],[0 1.2e4] * p(1) + p(2),'k');
[r,p]=corr(pds(~isnan(pds)),pdf(~isnan(pds)));
text(.65,.45,sprintf('r = %1.3f',r),'fontsize',8,'fontweight','bold','fontangle','normal');
ylabel({'from behavior'})
xlabel({'from stimulus'})
axis square;
axis tight;
box off;

tstr{1} = 'Variance explained';
tstr{2} = 'Comparison of embeddings';
tstr{3} = 'Embedding distances';

% -------------- LETTERS
hp=.05;
hy=1.15;
deffont=8;
for j = [1:length(hs)]
	hp0=hp;
	hy0 = hy;
	if j==2
		hp0=.11;
	end
	
	hpos = hs{j}.Position;
	lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
	axes('position', lpos);
	text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
	axis([0 1 0 1]);
	axis off;
	
	axes('position', [lpos(1)+.02 lpos(2)-.005 lpos(3) lpos(4)]);
	text(0,0, tstr{j},'fontangle','normal','fontsize',10);
	axis([0 1 0 1]);
	axis off;
end

%%
print(fullfile(matroot,'suppStimEmbed.pdf'),'-dpdf');







