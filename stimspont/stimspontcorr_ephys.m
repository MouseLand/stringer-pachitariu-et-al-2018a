function stimspontcorr_ephys(matroot)

%load('/media/carsen/DATA2/neuropix_stimspont.mat');
load('/home/carsen/dm11/data/Share/Carsen/neuropix_stimspont_v2.mat');

%%
igood = [];
for k = 1:length(rec)
	t0 = min(rec(k).gtRes);
	t1 = max(rec(k).gtRes);
	spont = zeros(max(rec(k).gtClu),floor((t1-t0)*100),'uint32');
	sbins = 0:0.01:t1-t0;
	for j = 1:max(rec(k).gtClu)
		ns = histcounts(rec(k).gtRes(rec(k).gtClu==j)-t0, sbins);
		spont(j,:) = uint32(ns);
	end
	spont = spont';
	igood = [igood mean(spont,1)>0.01];
	rec(k).igood = mean(spont,1)>0.01;
	spont = spont(:,mean(spont,1)>0.01);
	NN  = size(spont,2);
	
	istim   = isample(isample<=700);
	resp0   = squeeze(mean(double(rec(k).S(55+[1:50], isample<=700, rec(k).igood)),1));
	
	A = compute_means(istim, resp0, 2, 1);
	% only keep neurons that have signal variance
	vexp = diag(corr(A(:,:,1),A(:,:,2)));
	isig = vexp > 0.025;
	
	rec(k).igood = find(rec(k).igood);
	rec(k).igood = rec(k).igood(isig);
	spont = spont(:,isig);
	
	rec(k).spont = spont;
end
%%
igood = logical(igood(:));
wineu = ineu(igood);

tbins = [5 20 50 120 240 500 1000];
tstart = [55 49 49 49 49 49 49];

grps = {[1],[2],[3 4 5 6], [7], [8 9 10 11], [12 13 14]};

%grps = {[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14]};
%%
%rv= NaN*zeros(2,14,3);
clf;
%ac=[];
%acs=[];
for tb = [1:length(tbins)]
	clf;
	tbin = tbins(tb);
	for k = 1:length(rec)
		nt      = size(rec(k).spont,1);
		spont_bin   = bin2d(double(rec(k).spont), tbin, 1);
		spont_bin00 = spont_bin;
		tbin0 = 1;
		spont_bin0   = bin2d(double(rec(k).spont), tbin0, 1);
		
		if tbin > 50
			spont_bin00 = squeeze(mean(reshape(...
				double(rec(k).spont(1:floor(nt/50)*50,:)), 50, floor(nt/50), []),1));
		end
		
		
		istim   = isample(isample<=700);
		resp0   = squeeze(mean(double(rec(k).S(tstart(tb)+[1:min(tbin,50)], isample<=700, rec(k).igood)),1));
		
		mu      = mean(spont_bin00,1);
		sd      = std(spont_bin00,1,1)+ 1e-6;
		
		resp0   = (resp0 - mu)./sd;
		mu      = mean(spont_bin,1);
		sd      = std(spont_bin,1,1)+ 1e-6;
		spont_bin = (spont_bin - mu)./sd;
		spont_bin0 = (spont_bin0 - mu)./sd;
		
		A = compute_means(istim, resp0, 2, 1);
		% only keep neurons that have signal variance
		vexp = diag(corr(A(:,:,1),A(:,:,2)));
		isig = vexp > -inf;
		disp(mean(vexp));
		resp0 = resp0(:,isig,:);
		A = A(:,isig,:);
		spont_bin = spont_bin(:,isig);
		spont_bin0 = spont_bin0(:,isig);
		sineu = wineu(isig);
		
		stim_corr0 = corr(A(:,:,1), A(:,:,2));
		stim_corr0 = (stim_corr0 + stim_corr0')/2;
		spont_corr0 = corr(spont_bin);
		[u,s] = eig(spont_corr0);
		ac(:,:,k,tb) = acf(spont_bin0 * u(:,1:10), 1000);
		
		%if k==10
		Aex{tb}{k} = A;
		spontex{tb}{k} = spont_bin;
		%end
		
		for sub = 1:2
			if sub==2
				[u, s] = eig(stim_corr0 + spont_corr0);
				u = real(u(:,1)); v = u(:,1);
				stim_corr0 = stim_corr0 - u * (u'*stim_corr0*v) * v';
				spont_corr0 = spont_corr0 - u * (u'*spont_corr0*v) * v';
				acs(:,k,tb) = acf(spont_bin0 * u(:,1), 1000);
			end
			stim_corr = stim_corr0 - diag(NaN*diag(stim_corr0));
			spont_corr = spont_corr0 - diag(NaN*diag(spont_corr0));
			ninds = ~isnan(spont_corr(:)) & ~isnan(stim_corr(:));
			stim_corr = stim_corr(ninds);
			spont_corr = spont_corr(ninds);
			my_subplot(5,6,2*(k-1)+sub);
			%plot(mu,sd,'.');
			if sub==1
				plot(spont_corr,stim_corr,'r.','markersize',1);
			else
				
				plot(spont_corr,stim_corr,'.','markersize',1);
			end
			%xlabel('spont corr');
			%ylabel('stim corr');
			axis square;
			box off;
			
			[r,p] = corr(spont_corr,stim_corr); %,'type','Spearman');
			
			rv(sub,k,tb) = r;
			if k>1
				title(sprintf('r=%2.2f', r));
			else
				title(sprintf('no subtraction of shared\nr=%2.2f', r));
			end
			%title(sprintf('%2.2f %2.2f %2.2f %2.2f', r,mean(stim_corr),mean(spont_corr), mean(vexp)));
			drawnow;
		end
	end
end
%%

matroot = '/media/carsen/DATA2/grive/10krecordings/spontResults';
save(fullfile(matroot, 'stimspont_ephys.mat'),'rv','Aex','spontex','ac','acs','tbins');


%%
load(fullfile(matroot,'stimspont_ephys.mat'));

%%
close all;
default_figure([1 1 8 4]);

%%
clf;

signrank(rv(1,:,1), rv(2,:,1),'tail','right')


dex = 10;

clf;
clear hs;
xh=.55;
yh = xh;
i=0;
tb = 1;
stim_corr0 = corr(Aex{tb}{dex}(:,:,1), Aex{tb}{dex}(:,:,2));
stim_corr0 = (stim_corr0 + stim_corr0')/2;
spont_corr0 = corr(spontex{tb}{dex});

[u, s] = eig(stim_corr0 + spont_corr0);
u = real(u(:,1)); v = u(:,1);
u = u * -1;
[~,isort]=sort(u);

cm=[1 .2 0; .2 1 0; .2 0 1];
clear hs;
for sub = 1:2
	
	if sub==2
		stim_corr0 = stim_corr0 - u * (u'*stim_corr0*v) * v';
		spont_corr0 = spont_corr0 - u * (u'*spont_corr0*v) * v';
		
	end
	
	i=i+1;
	hs{i}=my_subplot(2,4,1+4*(sub-1),[xh yh]);
	hs{i}.Position(1) = hs{i}.Position(1)-.01;
	imagesc(spont_corr0(isort,isort)-diag(diag(spont_corr0(isort,isort))),[-.5 .5])
	xlabel('neurons')
	ylabel('neurons')
	title('spontaneous correlations','fontweight','normal')
	axis square;
	box off;
	if sub==2
		text(-.3,1.45,'1D Shared component subtracted','fontsize',10,'FontAngle','italic');
	end
	
	h=colorbar;
	h.Position(1)=	h.Position(1) +.06;
	
	i=i+1;
	hs{i}=my_subplot(2,4,2+4*(sub-1),[xh yh]);
	hs{i}.Position(1) = hs{i}.Position(1)-.02;
	imagesc(stim_corr0(isort,isort)-diag(diag(stim_corr0(isort,isort))),[-.25 .25])
	xlabel('neurons')
	ylabel('neurons')
	title('stimulus correlations','fontweight','normal')
	axis square;
	box off;
	h=colorbar;
	h.Position(1)=	h.Position(1) +.06;
	
	stim_corr = stim_corr0 - diag(NaN*diag(stim_corr0));
	spont_corr = spont_corr0 - diag(NaN*diag(spont_corr0));
	ninds = ~isnan(spont_corr(:)) & ~isnan(stim_corr(:));
	stim_corr = stim_corr(ninds);
	spont_corr= spont_corr(ninds);
	
	i=i+1;
	hs{i}=my_subplot(2,4,3+4*(sub-1),[xh yh]);
	hs{i}.Position(1) = hs{i}.Position(1)-.03;
	%plot(mu,sd,'.');
	plot(spont_corr,stim_corr,'.','color',cm(tb,:),'markersize',1);
	xlabel('spontaneous correlations');
	ylabel('stimulus correlations');
	axis square;
	[r,p] = corr(spont_corr,stim_corr);
	box off;
	b = polyfit(spont_corr, stim_corr, 1);
	hold all;
	plot([-.2 .61], b(2) + [-.2 .61]*b(1), 'k-')
	axis tight;
	ylim([-.25 .4]);
	text(.5, 1,sprintf('r=%2.2f', r),'fontweight','normal','fontsize',8);
	
end


i=i+1;
hs{i}=my_subplot(1,4,4);
irec=1:14;
tbs={'50 ms bins','200 ms','500 ms'};
for tb=1:3
	plot(rv(1,irec,tb),rv(2,irec,tb),'.','color', cm(tb,:), 'MarkerSize',8);
	hold all;
	plot(mean(rv(1,irec,tb),2),mean(rv(2,irec,tb),2),'x','color', cm(tb,:), 'MarkerSize',13,'MarkerFaceColor',cm(tb,:));
	text(.05,1.1-tb*.1,tbs{tb},'color',cm(tb,:),'fontsize',8);
end
plot([-.1 .5],[-.1 .5],'k')
axis square;
axis tight;
xlabel('before subtraction');
ylabel('after subtraction');
title({'correlation between','spont and stim correlations'},'FontSize',8,'FontWeight','normal');
grid on;
box off;
colormap(redblue)

hp=.05;
hy=1.24;
deffont=8;
for j = [1:length(hs)]
	hp0=hp; hy0=hy;
	if j==7
		hp0 =0.05;
		hy0 = .88;
	else
		hp0=hp;
		hy0 = hy;
		
	end
	hpos = hs{j}.Position;
	lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
	axes('position', lpos);
	t=text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','FontName', 'Helvetica');
	axis([0 1 0 1]);
	axis off;
	
end

%%
print(fullfile(matroot,'suppSSephys.pdf'),'-dpdf');


%%
close all;
default_figure([1 1 8 4]);

%%
clf;
dex=10;
k=0;
twindow = 20000 + [1 100000];
for tb = [1 4]
	k=k+1;
	spont_corr0 = corr(spontex{tb}{dex});
	
	[u, s] = eig(spont_corr0);
	u = real(u(:,1)); v = u(:,1);
	u = u * sign(skewness(u));
	[~,isort]=sort(u);
	
	my_subplot(2,1,k);
	imagesc(zscore(my_conv2(spontex{tb}{dex}(round(twindow(1)/tbins(tb)):round(twindow(2)/tbins(tb)), ...
		isort)', 1,1),1,2),[-0.5 2]);
	cm=colormap('gray');
	colormap(flipud(cm));
	title(tbins(tb)*10)
end


%%
clf;
clear hs;
cm = colormap('gray');
cm = cm([1 18 30 40 50],:);
colormap('redblue');
tbs={'50 ms','200 ms', '500 ms','1.2 s','2.4 s','5 s','10 s'};
dex=10;
for tb=1:length(tbins)
	stim_corr0 = corr(Aex{tb}{dex}(:,:,1), Aex{tb}{dex}(:,:,2));
	stim_corr0 = (stim_corr0 + stim_corr0')/2;
	spont_corr0 = corr(spontex{tb}{dex});
	
	[u, s] = eig(spont_corr0);
	u = real(u(:,1)); v = u(:,1);
	u = u * sign(skewness(u));
	if tb<3
		u = u*-1;
	end
	[~,isort]=sort(u);
	
	hp=my_subplot(2,7,tb,[.75 .75]);
	hp.Position(1) = hp.Position(1) + (5-tb)*.005 + .01;
	hp.Position(2) = hp.Position(2) + .01;
	if tb==1
		hs{1}=hp;
	end
	
	imagesc(spont_corr0(isort,isort)-diag(diag(spont_corr0(isort,isort))),[-.5 .5])
	title(tbs{tb},'fontweight','normal')
	if tb==1
		title(sprintf('time bin size = %s',tbs{tb}),'fontweight','normal')
	end
	if tb==1
		ylabel('neurons');
		xlabel('neurons');
		text(-.2,1.45,'spontaneous correlations','fontweight','normal','fontsize',10,'FontAngle','italic');
	end
	axis square;
	box off;
	if tb>1
		axis off;
	end
	
	
	hp=my_subplot(2,7,tb+7,[.7 .8]);
	hp.Position(1) = hp.Position(1) + (5-tb)*.005+.01;
	hp.Position(2) = hp.Position(2) + .01;
	if tb==1
		hs{2}=hp;
	end
	for n = 5:-1:1
		semilogx([0:size(ac,1)-1]*.01, squeeze(mean(ac(:,n,:,tb),3)),'color',cm(n,:));
		hold all;
		if tb==1
			text(.68,1-n*.1,sprintf('PC %d',n),'color',cm(n,:),'fontsize',8);
		end
	end
	%axis square;
	box off;
	ylim([-.05 .75]);
	xlim([0.01 10]);
	set(gca,'xtick',10.^[-2:1],'xticklabel',{'0.01','0.1','1','10'});
	title(tbs{tb},'fontweight','normal')
	if tb==1
		title(sprintf('time bin size = %s',tbs{tb}),'fontweight','normal')
	end
	
	if tb==1
		xlabel('time (s)');
		ylabel('autocorrelation');
		text(-.2,1.23,'autocorrelation of PC''s','fontweight','normal','fontsize',10,'FontAngle','italic');
	end
	grid on;
	grid minor;
	grid minor;
end

% k=0;
% ts = {'10 ms', '10 s'};
% for d =[2 1000]
% 	k=k+1;
% 	hp=my_subplot(3,4,9+k,[.6 .6]);
% 	hp.Position(2) = hp.Position(2) - .02;
% 	if k==1
% 		hs{3}=hp;
% 	end
% 	semilogx(tbins*.01,squeeze(mean(ac(d,n,:,:),3)),'color',cm(n,:));
% 	%if k==1
% 		xlabel('time bin size (s)')
% 		ylabel('autocorrelation');
% 	%end
% 	t(k)=title(sprintf('PC %d autocorr @ %s',n,ts{k}),'color',cm(n,:),'FontWeight','normal','fontsize',10,'FontAngle','italic');
% 	%t(k).Position(3)=0.6;
% 	axis square;
% 	box off;
% 	grid on;
% 	grid minor;
% 	grid minor;
% 	set(gca,'xtick',10.^[-2:1],'xticklabel',{'0.01','0.1','1','10'});
% end


hp=.04;
hy=1.05;
deffont=8;
for j = [1:length(hs)]
	hp0=hp; hy0=hy;
	if j>1
		hy0=1.22;
	end
	if j==3
		hp0=.06;
		
	end
	
	hpos = hs{j}.Position;
	lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
	axes('position', lpos);
	t=text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','FontName', 'Helvetica');
	axis([0 1 0 1]);
	axis off;
	
end


%%
print(fullfile(matroot,'suppCorrEphys.pdf'),'-dpdf');


%%
clf;
k=0;
for tb=[1 5]
	k=k+1;
	for n = 1:5
		my_subplot(2,5,5*(k-1)+n);
		semilogx([1:200]*.25, squeeze(ac(1:200,n,2,tb)),'r');
		hold all;
		semilogx([1:200]*.25, squeeze(acs(:,2,tb)),'b')
		axis square;
		axis tight;
		box off;
		legend(sprintf('PC %d @ %s',n,tbs{tb}),'shared');
		legend boxoff;
		xlabel('time (s)');
		set(gca,'xtick',[0.1 0.5 1  10]);
		grid on;
		grid minor;
		grid minor;
		if n==1
			title(tbs{tb})
		end
	end
	
end



%%

clf;
semilogx([0:size(ac,1)-1]*.01, squeeze(mean(ac(:,1:5,:,1),3)));
hold all;
semilogx([0:size(ac,1)-1]*.01, squeeze(mean(acs(:,:,tb),2)),'k')




















