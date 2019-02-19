function suppFastRasters(ephysroot, matroot)

load(fullfile(ephysroot, 'probeLocations.mat'));

areaLabels = {'FrCtx','FrMoCtx','SomMoCtx','SSCtx','V1','V2','RSP',...
    'CP','LS','LH','HPF','TH','SC','MB'};

% map from KS2 result order to Nick's probe location order
premap = [1 2 7 8 3 4 5 6];

mstr = {'Krebs','Waksman','Robbins'};
tstart = [3811 3633 3323]; % start of spontaneous activity

imouse = 2;
probeLoc = probeLocations(imouse).probe;
probeLoc = probeLoc(premap); % remap probes in same order as spks
mouse_name = mstr{imouse};

% load files
beh = load(fullfile(ephysroot, sprintf('%s_face_proc.mat',mouse_name)));
motSVD = beh.motionSVD;
tVid = beh.times; % times of movie frames in spike reference frame

load(fullfile(ephysroot, sprintf('spks%s_Feb18.mat',mouse_name)));

%% extract spikes
stall = zeros(5e3,1200000,'uint8');
ij = 0;
maxt=0;
Wh = [];
iprobe=[];
brainLoc=[];
srate = 500; % sampling rate in Hz (2 ms bins)
for k = 1:numel(spks)
	%%
	clu = spks(k).clu;
	st  = spks(k).st;
	st = round(st*srate);
	
	S = sparse(st, clu, ones(1, numel(clu)));
	S = uint8(full(S))';
	S = S(sort(unique(clu)),:);
	goodcells = mean(S,2)*srate > 0;
	S = S(goodcells,:);
	stall(ij+[1:size(S,1)],1:size(S,2)) = S;
	ij = ij + size(S,1);
	maxt = max(maxt, size(S,2));
	
	% height of spikes
	whp = spks(k).Wheights(sort(unique(clu)));
	whp = whp(goodcells);
	
	Wh = [Wh; whp];
	
	% where is the probe
	loc = zeros(numel(whp),1);
	lowerBorder = probeLoc(k).borders.lowerBorder;
	upperBorder = probeLoc(k).borders.upperBorder;
	acronym     = probeLoc(k).borders.acronym;
	%
	for j = 1:numel(acronym)
		whichArea = find(strcmp(areaLabels, acronym{j}));
		loc(whp >= lowerBorder(j) & whp < upperBorder(j)) = whichArea;
	end
	brainLoc = [brainLoc; loc];
	iprobe=[iprobe; k * ones(size(S,1),1)];
end

%
stall = stall(1:ij, 1:maxt);
%
tspont = tstart(imouse)*srate : min(floor(tVid(end)*srate), size(stall,2)-4);
stall = stall(:,tspont);

tspont = tspont / srate;

size(stall)

%% behavior
tBeh  = tspont;
x = interp1(tVid, motSVD, tBeh);

stall = stall';

%% 1st PC and manifold embedding sorting
nC =30;
ops.nCall = [nC 100]; % number of clusters
ops.iPC = 1:200; % PCs to use
ops.useGPU = 1; % whether to use GPU
ops.upsamp = 100; % upsampling factor for the embedding position
ops.sigUp = 1; % stddev for upsampling
S = bin2d(stall', 100, 2);
S = zscore(S,1,2) ./ sqrt(size(S,2));

[isort, ~, Sm] = mapTmap(S, ops);

isort = gather_try(isort);

%%
sr = .3 * srate;
spont_bin = bin2d(single(stall), sr, 1);

%%
grps = {[1 2], [3 4], [5 6], [7],[8],[9 10 13 14], [11],[12]};
tgrps = {'frontal','sensorimotor','visual','retrosplenial','striatum','midbrain','hippocampus','thalamus'};
cg = colormap('hsv');
cg = cg(1:8:end,:);
cg = cg(1:8,:);
cg = max(0,cg-.1);
rng(1);
cg = cg(randperm(8),:);
%cg = [0 0 0; cg];
ngrp = zeros(numel(grps),1);
ng = numel(grps);
iarea = zeros(length(isort),1);
for j = 1:numel(grps)
	for k = grps{j}
		iarea(brainLoc==k) = j;
	end
end
isort = isort(iarea(isort)>0);

i=0;
%%
close all;
default_figure([1 1 7.25 7.25]);
%%
tgrps = {'frontal',{'sensori-','motor'},'visual',{'retro-','splenial'},'striatum','midbrain','hippocampus','thalamus'};
i1=5e4/2000;
%i1=ceil(rand*100)
i1=60; % 12, ,4, 105
i2 = 250; % 60 , 56, 20, 200, 302, 430, 462
i2 = i2(end);
clf;
cm=colormap('gray');
colormap(flipud(cm));
twindow = [(i1-8)*2000/sr (i2+8)*2000/sr];%[1 size(spont_bin,1)];
%tw = 2.2e4 + [1 4000];
tw = {i1*2000 + [1 2000], i2*2000 + [1 2000]};

clear hs;
hp=my_subplot(3,1,1,[.75 .45]);
hp.Position(1) = hp.Position(1) - .035;
imagesc(zscore(my_conv2(spont_bin(:, isort)', [3 2], [1 2]),1,2),[-0.5 2]);
hold all;
for tt=1:length(tw)
	patch([tw{tt}(1) tw{tt}(1) tw{tt}(2) tw{tt}(2)]/sr,[-10 length(isort)+10 length(isort)+10 -10],'k',...
		'facecolor','none','edgecolor',[1 .7 .7],'linewidth',1);
end
axis tight;
%imagesc(zscore(my_conv2(single(stall(tw(1) : tw(2), isort(iarea(isort)==j))'),[1 3],[1 2]),1,2),[0 0.05])
axis off;
xlim([twindow(1) twindow(2)]);

axes('position',[hp.Position(1) hp.Position(2)+hp.Position(4)+.005 hp.Position(3) .05]);
xt = x(:, 1:3);
xt =xt*-1;
xt=xt-mean(xt,1);
for j=1:3
	plot(my_conv2(xt(:,j),20,1),'color',[.3 j*0.3 (4-j)*.3]);
	hold all;
end
text(0,1,'face motion PCs','color',[.3 .3 .9]);
%plot(x(tw{tt}(1) : tw{tt}(2), 2) *-1,'r');
%plot(,'b');
axis off;
text(-0.03, 1.15,'A','fontweight','bold','fontsize',10);
xlim([twindow(1)*sr twindow(2)*sr]);
ylim([-100 1500]);

axes('position',[hp.Position(1)-.002 hp.Position(2)-.0002 hp.Position(3) hp.Position(4)]);
plot([0 10/(sr/srate) ],[0 0 ],'k');
hold all;
plot([0 0], [0 500],'k');
xlim([0 twindow(2)-twindow(1)]);
ylim([0 size(spont_bin,2)]);
text(-0.01,0,'10 s');
ht=text(-.02,0,'500 neurons');
set(ht,'rotation',90);
axis off;

for tt=1:2
hp2=my_subplot(1,2,tt,[.8 .68]);
hp2.Position(2) = hp2.Position(2) - .19;
if tt==1
	hpmarge = hp.Position(1) - hp2.Position(1);
	hpsec = hp.Position(3) - hp2.Position(3) + 0.02;%hpmarge/2;
else
	hpmarge = (hp2.Position(1)+hp2.Position(3)) - (hp.Position(1)+hp.Position(3));
	hpsec = hp.Position(3) - hp2.Position(3)- .01;% hpmarge/2;
end

axis off;
dy=0.04;
axes('position',[hp2.Position(1) hp2.Position(2)+hp2.Position(4)+.02 hp2.Position(3) hp.Position(2)-hp2.Position(2)-hp2.Position(4)-.02]);
hold all;
ttot = diff(twindow);
tsub = ttot*hpmarge;
if tt==1
	plot([tw{tt}(1)/sr twindow(1)-tsub], [10 0],'color',[1 .7 .7]);
	plot([tw{tt}(2)/sr twindow(1)+hpsec*diff(twindow)], [10 0],'color',[1 .7 .7]);
else
	plot([tw{tt}(2)/sr twindow(2)+tsub], [10 0],'color',[1 .7 .7]);
	plot([tw{tt}(2)/sr twindow(2)-hpsec*diff(twindow)], [10 0],'color',[1 .7 .7]);
end
%	plot([mean(tw{tt})/sr twindow(2)-diff(twindow)*.15], [10 0],'color',[1 .7 .7]);
%end
axis tight;
if tt==1
	xlim([twindow(1)-tsub twindow(1) + hpsec * diff(twindow)]);
else
	xlim([twindow(2) - hpsec * diff(twindow) twindow(2) + tsub]);
end
axis off;
end

NN = size(stall,2);
isortneu = [];
ij =[];
sts{1}=[];
sts{2}=[];
for j = [3 4 2 1 7 8 6]%unique(iarea(iarea>0))'
	ineu = isort(iarea(isort)==j);
	ineu = ineu(end:-1:1);
	%ineu = ineu(randperm(length(ineu)));
	%[~,ix] = sort(Wh(ineu));
	%ineu = ineu(ix);
	for tt = 1:2
		sts{tt} = cat(2,sts{tt}, my_conv2(single(stall(tw{tt}(1) : tw{tt}(2), ineu)), [3 1],[1 2]));
		sts{tt} = cat(2,sts{tt}, zeros(tw{tt}(2)-tw{tt}(1)+1, 30, 'single'));
	end
	ineu = [ineu; NN * ones(30,1)];
	isortneu = [isortneu; ineu(:)];
	ij = [ij; ones(length(ineu),1)*j];
	
end
ts = {'B','C'};
for tt=1:2
	hp=my_subplot(1,2,tt,[.8 .68]);
	hp.Position(2) = hp.Position(2) - .19;
	imagesc(zscore(sts{tt}',1,2),[0 0.2])
	text(-.1, 1.05,ts{tt},'fontweight','bold','fontsize',10);
	jdiff = find(abs(diff(ij))>0);
	hold all;
	k=0;
	jdiff = [0;jdiff;length(isortneu)];
	for k = 2:length(jdiff)
		%plot([1 tw{tt}(2)-tw{tt}(1)], [jdiff(k) jdiff(k)], 'k');
		ar = ij(jdiff(k));
		ht=text(-200,(jdiff(k)-jdiff(k-1))/2 + jdiff(k-1),tgrps{ar},'color',cg(ar,:),'fontsize',7,'fontweight','bold',...
			'HorizontalAlignment','center','units','data');
		set(ht,'rotation',90);
	end
	hold all;
	plot([-20 250 ],size(sts{tt},2)-15+[0 0 ],'k');
	plot([0 0]-20, size(sts{tt},2)-15+[0 -115],'k');
	text(0,0.004,'500 ms');
	ht=text(-.04,0,'100 neurons');
	set(ht,'rotation',90);
	axis tight;
	axis off;
	
	axes('position',[hp.Position(1) hp.Position(2)+hp.Position(4) hp.Position(3) .04]);
	xt = x(tw{tt}(1) : tw{tt}(2), 1:3);
	xt =xt*-1;
	xt=xt-mean(xt,1);
	for j=1:3
		plot(xt(:,j),'color',[.3 j*0.3 (4-j)*.3]);
		hold all;
	end
	text(0.4,.55*(tt-1)+.8,'face motion PCs','color',[.3 .3 .9]);
	axis tight;
	axis off;
	xlim([-20 2000]);
	ylim([-500 1400]);
	
end


%%
%%
print(fullfile(matroot,'suppEphysEx2.pdf'),'-dpdf');