function fig3new(matroot)

%load('exampleFrame.mat');
load(fullfile(matroot, 'ephysResults_KS2.mat'));
%%
% diagram with probes
close all;
default_figure([1 1 4.75 4.]);
%%
i=0;
clf;
dex = 2;
cdd=colormap('cool');
cdd = cdd(35:5:end,:);
cm =colormap('gray');

xh = .6;
yh = xh;
clear hs;
i=i+1;
hs{i}=my_subplot(3,4,1,[.8 .8]);
%hs{i}.Position(1) = hs{i}.Position(1) - .02;
%hs{i}.Position(2) = hs{i}.Position(2) - .06;
im=imread('allenprobes.png');
imagesc(im);
axis off;
axis image;
pos = hs{i}.Position;
i=i+1;
hs{i}=my_subplot(3,4,2,[xh yh]);
hs{i}.Position(1) = hs{i}.Position(1) - .03;
im=imread('image.png');
imagesc(im);
axis off;
axis image;
text(1,1,'calbindin','color',[0 1 0],'fontweight','bold','fontangle','normal','HorizontalAlignment','right','fontsize',6);
text(1,.85,'dil','color',[1 0 0],'fontweight','bold','fontangle','normal','HorizontalAlignment','right','fontsize',6);
%title('percent of neurons');

%title('where are the probes (cross-section)');
grps = {[1 2], [3 4], [5 6], [7],[8],[9 10 13 14], [11],[12]};
tgrps = {'frontal','sensorimotor','visual','retrosplenial','striatum','midbrain','hippocampus','thalamus'};
cg = colormap('hsv');
cg = cg(1:8:end,:);
cg = cg(1:8,:);
cg = max(0,cg-.1);
rng(1);
cg = cg(randperm(8),:);
cg = [0 0 0; cg];
ngrp = zeros(numel(grps),1);
ng = numel(grps);
clear iarea;
for d = 1:3
    iarea{d} = zeros(numel(results.isort{d}),1);
end

for j = 1:numel(grps)
    for k = grps{j}
        ngrp(j) = ngrp(j) + results.allLoc(k);
        for d = 1:3
            iarea{d}(results.brainLoc{d}==k) = j;
        end
    end
   
end
set(gca,'fontsize',6);

% i=i+1;
% hs{i}=my_subplot(3,4,2,[xh yh]);
ng = length(grps);
sp = 0.03;

% for j = 1:length(ngrp)
% 	hold all;
% 	bar(j,ngrp(j),'facecolor',cg(j+1,:));
% 	ht=text(xg(j)-.01,-.02,tgrps{j},'HorizontalAlignment','right','fontsize',7,'fontangle','normal','color',cg(j+1,:));
% 	ht.Rotation = 45;
% 	%bar(ngrp,'color',cg);
% end
% ylabel('# of neurons');
% box off;
% set(gca,'xtick',[]);
% %axis square;
% axis tight;


i=i+1;
hs{i}=my_subplot(3,4,3,[xh yh]);
cin = zeros(ng,1);
cout = zeros(ng,1);
for d = 1:3
    cc = corr(results.spks{d}');
    c2 = triu(cc);
    itri = find(abs(c2)~=0);
    cc2 = cc;
    cc2(itri) = NaN;
    cc = cc-diag(NaN*diag(cc));
    for j = 1:ng
        if d==1
            cain{j} = [];
            caout{j} = [];
        end
        csub = (cc2(iarea{d}==j,iarea{d}==j));
        cain{j} = [cain{j}; csub(~isnan(csub))];
        csub = (cc(iarea{d}==j,iarea{d}~=j));
        caout{j} = [caout{j}; csub(~isnan(csub))];
    end
end
hold all;
plot([.0 .11],[.0 .11],'k');
for j = 1:ng
    plot(mean(cain{j}), mean(caout{j}),...
        '.','color',cg(j+1,:),'markersize',10);
	
end
set(gca,'xtick',[0:.05:.15],'ytick',[0:.05:.15]);
%axis([ 0 .12]);
box off;
%axis([.5 8.5 -.15 .45]);
axis square;
axis tight;
ylabel({'out-of-area corr'});
xlabel({'within-area corr'});
set(gca,'fontsize',6);

% ---------- CORR WITH FIRST PC ---------- %
i=i+1;
hs{i}=my_subplot(3,4,4,[xh yh]);
for d = 1:3
    cc = results.corr_pc1{d};%corr(results.spks{d}',results.pc1{d}(:));
    for j = 1:ng
        if d==1
            ca1{j}=[];
        end
        ca1{j} = [ca1{j}; cc(iarea{d}==j)];
    end
end
box off;
set(gca,'xtick',[1:8],'xticklabel',{});
xtickangle(45);
ylim([-.75 .75]);
ylabel({'corr with face PC1'});
xg = linspace(sp,1-sp,ng);
perms = [4 3 1 2 5 8 6 7];
xg = xg(perms);
for j = 1:ng
	hold all;
	errorbar(perms(j), mean(ca1{j}), std(ca1{j}),'.', 'markerfacecolor', cg(j+1,:),...
		'markersize',10,'color', cg(j+1,:))
	ht=text(xg(j)-.01,-.02,tgrps{j},'HorizontalAlignment','right','fontsize',6,'fontangle','normal','color',cg(j+1,:));
	ht.Rotation = 45;
end
%axis square;
%axis([.5 8.5 -.38 .48]);
%axis tight;
xlim([.5 8.3]);
grid on;
set(gca,'fontsize',6);
axis square;

% --------------- RASTER -------------%
i=i+1;
hs{i}=my_subplot(3,2,3,[0.9 .78]);
hs{i}.Position(1) = hs{i}.Position(1) + .0;
hs{i}.Position(2) = hs{i}.Position(2) + .01;
trange = [1:800];
NN=size(results.spks{dex},1);
%imagesc(zscore(ypred(gather(isort), :),1,2),[0 .3]);
imagesc(my_conv2(zscore(results.ytest{dex}(results.isort{dex}(end:-1:1),trange),1,2),2,1),[0.0 .5]);
cm = colormap('gray');
colormap(flipud(cm));
hold all;
plot([0 0]-5,NN+[-500 0],'k');
ht=text(-.05,0,'500 neurons','fontangle','normal','fontsize',6);
ht.Rotation=90;
plot([-5 60/1.2],(NN+10)*[1 1],'k');
text(0.0,-.03,'1 min','fontsize',6,'fontangle','normal');
xlim([-2 numel(trange)]);
axis tight;
axis off;
title('       Neural activity (test set) sorted by 1D embedding','fontweight','normal','fontsize',6);
set(gca,'fontsize',6);

hp=my_subplot(3,2,5,[.9 .78]);
hp.Position(2) = hp.Position(2) - .0;
hp.Position(1) = hs{i}.Position(1);
NN=size(results.spks{dex},1);
imagesc(zscore(results.ypred{dex}(results.isort{dex}(end:-1:1),trange),1,2),[0.0 .5]);
hold all;
plot([0 0]-5,NN+[-500 0],'k');
ht=text(-.05,0,'500 neurons','fontangle','normal','fontsize',6);
ht.Rotation=90;
plot([-5 60/1.2],(NN+10)*[1 1],'k');
text(0.0,-.03,'1 min','fontsize',6,'fontangle','normal');
xlim([-2 numel(trange)]);
axis tight;
axis off;
title('Neural activity prediction (test set) from face motion','fontweight','normal','fontsize',6);
set(gca,'fontsize',6);

% neurons along continuum
hp=my_subplot(3,4,7,[xh*1.2 yh]);
hp.Position(2) = hs{i}.Position(2);
hp.Position(4) = hs{i}.Position(4);
hp.Position(1) = hp.Position(1)-.05;
isort = results.isort{dex};
wh = 4000-results.Wh{dex}(isort);
hold all;
for j = 0:numel(grps)
    wp=wh(iarea{dex}(isort)==j)+randn(sum(iarea{dex}(isort)==j),1);
    ip = find(iarea{dex}(isort)==j);
    plot(wp,ip,'.','color',cg(j+1,:),'markersize',2);
end 
xlabel({'depth (mm)'});
ylabel('embedding position');
set(gca,'ytick',[],'YAxisLocation','right','xtick',[0:2000:4000],'xticklabel',{'0','2','4'});%'xaxisLocation','top',
box off;
axis tight;
xlim([0 4000]);
set(gca,'fontsize',6);



% -------------- BEHAVIOR ---- %
%idelay = find(results.tdelay==0.20);
idelay = 19;
i=i+1;
hs{i}=my_subplot(3,4,8,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) - .04;
hs{i}.Position(1) = hs{i}.Position(1) + .01;
bpred = [];
npred = [];
for j = 1:ng
    for d = 1:3
        if d==1
            resbeh{j}=[];
            totbeh{j}=[];
			resall{j}=[];
            totall{j}=[];
        end
        resbeh{j} = [resbeh{j}; results.resid_behavior{d}(iarea{d}==j,:,:)];
        totbeh{j} = [totbeh{j}; results.totvar_behavior{d}(iarea{d}==j,:,:)];
		resall{j} = [resall{j}; results.resid{d}(iarea{d}==j,:)];
        totall{j} = [totall{j}; results.totvar{d}(iarea{d}==j,:)];
    end
    %semilogx(nPC, 1 - nanmean(resall{j})./nanmean(totall{j}), 'color', cg(j+1,:));
    semilogx(results.ndims0, 100 * (1-nanmean(resbeh{j}(:,:,idelay))./nanmean(totbeh{j}(:,:,idelay))),...
        'color', cg(j+1,:));
	bpred(j) = 100 * (1-nanmean(resbeh{j}(:,end,idelay))./nanmean(totbeh{j}(:,end,idelay)));
	npred(j) = 100 * (1-nanmean(resall{j}(:,end-1))./nanmean(totall{j}(:,end-1)));
	fprintf('%s beh: %2.1f neur: %2.1f ratio: %2.1f \n',tgrps{j}, bpred(j), npred(j), bpred(j) / npred(j) * 100);
	
    hold all;
end
axis([1 32 0 37]);
set(gca,'xtick',2.^[0:2:10]);
ylabel({'% variance explained'});
xlabel('rank');
box off;
grid on;
grid minor;
grid minor;
axis square;
set(gca,'fontsize',6);

% pos = hs{i}.Position;
% axes('position',[pos(1)+pos(3)-.015 pos(2)+pos(4)*.5 pos(3)*.5 pos(4)*.5]);
% imagesc(fr(100:235,200:335),[0 150]);
% colormap(gca, cm);
% axis image;
% axis off;

avev = 1-nanmean(results.resid_behavior{1}(:,5,idelay))/nanmean(results.totvar_behavior{1}(:,5,idelay));

bb = 1-nanmean(resbeh{j}(:,:,idelay))./nanmean(totbeh{j}(:,:,idelay));


% ------------ TIMELAGS------------%
i=i+1;
hs{i}=my_subplot(3,4,11,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) - .02;
hs{i}.Position(1) = hs{i}.Position(1) + .02;
hold all;
for j = 1:ng
    plot(results.tdelay(1:end), 100* squeeze(1-nanmean(resbeh{j}(:,5,:))./nanmean(totbeh{j}(:,5,:))),...
        'color', cg(j+1,:),'linewidth',1);
    hold all;
end
box off;
ylabel({'% variance explained'});
xlabel('time from behav. (s)');
grid on;
grid minor;
grid minor;
axis square;
set(gca,'xtick',[-5 0 5],'xticklabel',{'-5','0 ','5'});
axis ([-5 5 0 36]);
set(gca,'fontsize',6);

i=i+1;
hs{i}=my_subplot(3,4,12,[xh yh]);
hs{i}.Position(2) = hs{i}.Position(2) - .02;
hs{i}.Position(1) = hs{i}.Position(1) + .01;
resid = [];
totvar=[];
for d = 1:3
    resid = cat(1,resid, results.resid_bshort{d});
    totvar = cat(1,totvar, results.totvar_bshort{d});
end
hold all;
for j = 1:ng
    for d = 1:3
        if d==1
            resbeh{j}=[];
            totbeh{j}=[];
        end
        resbeh{j} = [resbeh{j}; results.resid_bshort{d}(iarea{d}==j,:,:)];
        totbeh{j} = [totbeh{j}; results.totvar_bshort{d}(iarea{d}==j,:,:)];
	end
	plot(results.tdelay(1:end), 100* squeeze(1-nanmean(resbeh{j}(:,5,:))./nanmean(totbeh{j}(:,5,:))),...
        'color', cg(j+1,:),'linewidth',1);
    hold all;
end
% for d = [1:5]
%     plot(results.tdelay, squeeze(1-nanmean(resid(:,d,:))./...
%         nanmean(totvar(:,d,:)))',...
%         'color', cdd(d,:),'linewidth',0.5);
%     text(.8, d/8+.4, sprintf('%dD',results.ndims0(d)),...
%         'color', cdd(d,:),'fontangle','normal','horizontalalignment','right','fontsize',6);
% 	re(:,d) = squeeze(1-nanmean(resid(:,d,:))./...
%         nanmean(totvar(:,d,:)))';
% end
box off;
ylabel({'% variance explained'});
xlabel('time from behav. (s)');
grid on;
axis tight;
ylim([0 18]);
grid minor;
grid minor;
axis square;
set(gca,'xtick',[-5 0 5],'xticklabel',{'-5','0 ','5'});
set(gca,'fontsize',6);
%axis ([-5 5 -.005 .05]);
xlim([-5 5]);
%
resid = [];
totvar=[];
for d = 1:3
    resid = cat(1,resid, results.resid_behavior{d});
    totvar = cat(1,totvar, results.totvar_behavior{d});
end
[~,bestdelay]=max(1-squeeze(resid(:,5,:))./squeeze(totvar(:,5,:)),[],2);
md = mean(results.tdelay(bestdelay));
stdd=std(results.tdelay(bestdelay))/sqrt(numel(bestdelay)-1);
disp([md stdd]);

%
tstr{1} = 'probe locations';
tstr{2} = '';
tstr{3} = '';
tstr{4} = '';
tstr{5} = '';
tstr{6} ='from face motion';
tstr{7} ='time lag analysis';
tstr{8} =  '200 ms time bins';

% -------------- LETTERS
hp=.08;
hy=1.25;
deffont=6;
for j = [1:length(hs)]
	hy0 =hy;
	hp0=hp;
    if j==1 || j==2
        hp0 =.02;    
		hy0 = 1.05;
		if j==2
			hy0=1.2;
		end
	elseif j==3 || j==4
		hy0 = 1.2; 
	elseif j==5
		hp0=0.02;
		hy0 = 1.16;
    end
    hpos = hs{j}.Position;
    lpos = [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01];
    axes('position', lpos);
    text(0,0, char(64+j),'fontsize',deffont+4,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
    
    axes('position', [lpos(1)+.03 lpos(2)-.01 lpos(3) lpos(4)]);
    text(0,0, tstr{j},'fontangle','normal','fontsize',7);
    axis([0 1 0 1]);
    axis off;
    
end


%%
%
print(fullfile(matroot,'fig3new.pdf'),'-dpdf');


%% 1.2 second bin fwhm

tl = 1-squeeze(mean(resid(:,5,:),1))./squeeze(mean(totvar(:,5,:),1));

tli = interp1(results.tdelay,tl,[-8:.01:8]);

[~,ix1]=min(abs(max(tli)/2 - tli(1:800)));
[~,ix2]=min(abs(max(tli)/2 - tli(801:end)));

fwhm = (800-ix1+ix2)*.01

%% short timescales fwhm
resid = [];
totvar=[];
for d = 1:3
    resid = cat(1,resid, results.resid_bshort{d});
    totvar = cat(1,totvar, results.totvar_bshort{d});
end

for nd = 1:5
    tl = 1-squeeze(mean(resid(:,nd,:),1))./squeeze(mean(totvar(:,nd,:),1));

    tli = interp1(results.tdelay,tl,[-8:.01:8]);

    [~,ix1]=min(abs(max(tli)/2 - tli(1:800)));
    [~,ix2]=min(abs(max(tli)/2 - tli(801:end)));

    fwhm = (800-ix1+ix2)*.01
end



%%
drawnow;
pause(3);


%% --------- PEERS ---------%
clf;
nPC = 2.^[0:9];
expv_neurons = zeros(ng, 10);
for j = 1:ng
    for d = 1:3
        if d==1
            resall{j}=[];
            totall{j}=[];
        end
        resall{j} = [resall{j}; results.resid{d}(iarea{d}==j,:)];
        totall{j} = [totall{j}; results.totvar{d}(iarea{d}==j,:)];
    end
    semilogx(nPC, 1 - nanmean(resall{j})./nanmean(totall{j}), 'color', cg(j+1,:));
	expv_neurons(j,:) = (1 - nanmean(resall{j})./nanmean(totall{j}));
    hold all;
end
set(gca,'xtick',2.^[0:3:10]);
axis([1 512 0 1]);
box off;
ylabel({'variance explained','(test set)'});
xlabel('dimensions');
grid on;
grid minor;
grid minor;
axis square;
colors = cg(2:end,:);
save(fullfile(matroot, 'ephys_peers.mat'), 'expv_neurons','colors','tgrps')
