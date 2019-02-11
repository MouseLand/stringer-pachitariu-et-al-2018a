function supppcs(matroot)

load(fullfile(matroot,'corr1stpc.mat'));

%%
close all;
default_figure([10 10 6 6]);
%%
clf;
ndat = length(results.autocorr);
for d = 1:ndat
    my_subplot(3,3,d);
    for j = 1:8
        semilogx([1:400]/1.2,results.autocorr{d}(:,j),'color',(j-1)/9*[1 1 1]);
        hold all;
        text(.75,1.1-j*.08,sprintf('PC%d',j),'color',(j-1)/9*[1 1 1]);
    end
    set(gca,'xtick',[1/1.2 10.^[1:2]],'xticklabel',{'0','10','100'});
    box off;
    axis tight;
    ylim([-.1 1]);
    axis square;
    xlabel('time (s)');
    grid on;
    grid minor;
    grid minor;
end

%%
print(fullfile(matroot,'suppPCAC.pdf'),'-dpdf','-bestfit');
