function stimspontcorr_ephys(matroot)

load('/media/carsen/DATA2/neuropix_stimspont.mat');

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
	rec(k).spont = spont;
end
%%
igood = logical(igood(:));
wineu = ineu(igood);

tbins = [5 25 50];
tstart = [55 49 49];

grps = {[1],[2],[3 4 5 6], [7], [8 9 10 11], [12 13 14]};

%grps = {[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14]};
%%
rv= NaN*zeros(2,14,3);
clf;
for tb = 1:3
	clf;
	tbin = tbins(tb);
	for k = 1:length(rec)
		nt      = size(rec(k).spont,1);
		spont_bin   = squeeze(mean(reshape(double(rec(k).spont(1:floor(nt/tbin)*tbin,:)), tbin, floor(nt/tbin), []),1));
		istim   = isample(isample<=700);
		resp0   = squeeze(mean(double(rec(k).S(tstart(tb)+[1:tbin], isample<=700, rec(k).igood)),1));
		
		mu      = mean(spont_bin,1);
		sd      = std(spont_bin,1,1)+ 1e-6;
		
		resp0   = (resp0 - mu)./sd;
		spont_bin = (spont_bin - mu)./sd;
		
		A = compute_means(istim, resp0, 2, 1);
		% only keep neurons that have signal variance
		vexp = diag(corr(A(:,:,1),A(:,:,2)));
		isig = vexp > 0.025;
		disp(mean(vexp));
		resp0 = resp0(:,isig,:);
		A = A(:,isig,:);
		spont_bin = spont_bin(:,isig);
		sineu = wineu(isig);
		
		stim_corr0 = corr(A(:,:,1), A(:,:,2));
		stim_corr0 = (stim_corr0 + stim_corr0')/2;
		spont_corr0 = corr(spont_bin);
		
		if tb==1 && k==10
			Aex=A;
			spontex=spont_bin;
		end
		
		for sub = 1:2
			if sub==2
				[u, s] = eig(stim_corr0 + spont_corr0);
				u = real(u(:,1)); v = u(:,1);
				stim_corr0 = stim_corr0 - u * (u'*stim_corr0*v) * v';
				spont_corr0 = spont_corr0 - u * (u'*spont_corr0*v) * v';
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

save(fullfile(matroot, 'stimspont_ephys.mat'),'rv','Aex','spontex');


%%
load(fullfile(matroot,'stimspont_ephys.mat'));

close all;
default_figure([1 1 8 4]);

%%
clf;

signrank(rv(1,:,1), rv(2,:,1),'tail','right')

stim_corr0 = corr(Aex(:,:,1), Aex(:,:,2));
stim_corr0 = (stim_corr0 + stim_corr0')/2;
spont_corr0 = corr(spontex);

[u, s] = eig(stim_corr0 + spont_corr0);
u = real(u(:,1)); v = u(:,1);
[~,isort]=sort(u);

clf;

xh=.55;
yh = xh;
i=0;
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
	plot(spont_corr,stim_corr,'r.','markersize',1);
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

cm=[1 .2 0; .2 1 0; .2 0 1];

i=i+1;
hs{i}=my_subplot(1,4,4);
irec=1:14;
tbs={'50 ms bins','250 ms','500 ms'};
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






















