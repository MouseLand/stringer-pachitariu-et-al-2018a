function suppVarPos(matroot)

load(fullfile(matroot,'expv_neurons_pos.mat'));

%%

% compute distances between high and low variance neurons
for d = 1:numel(cellpos)
    cpos = cellpos{d};
    ev = expv_neurons{d}(:,6);
    eg{1} = ev<.03;
    eg{2} = ev>.1;
    disp(sum([eg{1} eg{2}]));
    NN = size(cpos,1);
    pd0 = zeros(NN,NN);
    dp{1} = 3; % depth distances
    dp{2} = [1:2]; % XYZ distances
    for k = 1:2
        for j = dp{k}
            pd0 = pd0 + (cpos(:,j) - cpos(:,j)').^2;
        end
        pd = sqrt(pd0);
        jn = [2 1];
        for j = 1:2
            din(j,k,d) = nanmean(nanmean(pd(eg{j},eg{j})));
            dout(j,k,d) = nanmean(nanmean(pd(eg{j},eg{jn(j)})));
        end
    end
    
    % depth distribution of variance explained
    depths = unique(cpos(:,3));
    depths = sort(depths);
    for j =1:numel(depths)
        dev(j,d) = nanmean(ev(cpos(:,3)==depths(j)));
    end
end

% dataset colors
load('cdat.mat');
%%
close all;
default_figure([15 1 4 4]);

%%
load('cdat.mat');
ndat = size(cdat,1);
clf;
i=0;
clear hs;
i=i+1;
hs{i}=my_subplot(2,2,1,[.65 .65]);
axis off;
hp=hs{i}.Position;
hp(2) = hp(2)-.03;
axes('position',hp);
hold all;
for d = 1:ndat
    plot(depths,dev(:,d),'color',cdat(d,:),'linewidth',.5);
end
%ylim([0.05 .3]);
xlabel('depth (um)');
ylabel('variance explained');
box off;
xlim([70 35*12])
axis tight;
camroll(-90)
text(-.24,1.55, 'depth distribution');
text(1.7,1.55, 'per cluster depth distance');
text(3.7,1.55, 'depth distance averaged');

i=i+1;
hs{i}=my_subplot(2,2,2,[.65 .65]);
hold all;
for d = 1:ndat
    plot([din(:,1,d) dout(:,1,d)],'color',cdat(d,:));
end
axis tight;
axis([.75 2.25 0 130])
box off;
axis square;
ylabel('pairwise distance (\mum)')
set(gca,'xtick',[1 2],'xticklabel',{'same variance','diff variance'});
%xtickangle(45);

i=i+1;
hs{i}=my_subplot(2,2,3,[.65 .65]);
axis off;
my_subplot(2,2,3,[.75 .75]);
% plot example plane
dex = 9;
iplane = 5;
hold all;
cpos = cellpos{dex}(cellpos{dex}(:,3)==iplane*35,1:2);
ev   = expv_neurons{dex}(cellpos{dex}(:,3)==iplane*35, 6);
cpos(isnan(ev),:) = [];
ev(isnan(ev))=[];

evsort = sort(ev);
bed =evsort(round(linspace(1,numel(ev),100)));
[~,~,bins]=histcounts(ev,bed);
NN  = size(cpos,1);
hcol = .66*bins/max(bins);
hsv = [.66-hcol(:) ones(NN,1) ones(NN,1)];
cm = hsv2rgb(hsv);
for j = 1:NN
    if ev(j) > 1e-2
        plot(cpos(j,1),cpos(j,2),'.',...
            'color','k','markersize',8*(log10(ev(j))+2));
    else
        plot(cpos(j,1),cpos(j,2),'x',...
            'color',.5*[1 1 1],'markersize',2,'linewidth',.5);
    end
end
axis tight;
axis off;
axis square;


i=i+1;
hs{i}=my_subplot(2,2,4,[.65 .65]);
hold all;
for d = 1:ndat
    plot([din(:,2,d) dout(:,2,d)],'color',cdat(d,:));
end
axis tight;
ylim([0 580]);
xlim([0.75 2.25]);
box off;
axis square;
ylabel('pairwise distance (\mum)')
set(gca,'xtick',[1 2],'xticklabel',{'same variance','diff variance'});


% -------------- LETTERS
hp=.13;
hy=1.25;
deffont=8;
for j = [1:length(hs)]
	hp0=hp;
	hy0=hy;
    hpos = hs{j}.Position;
    axes('position', [hpos(1)-hp0 hpos(2)+hpos(4)*hy0 .01 .01]);
    text(0,0, char(64+j),'fontsize',deffont+2,'fontweight','bold','fontangle','normal');
    axis([0 1 0 1]);
    axis off;
end

%
print(fullfile(matroot,'suppVarV1.pdf'),'-dpdf','-bestfit');