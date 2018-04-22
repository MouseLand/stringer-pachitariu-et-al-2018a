function layeredNet(axpos, nls, call, msize, ells)

axes('position',axpos);

bcent = max(nls)/2;
hold all;

for k = 1:numel(nls)-1
    nl1 = nls(k);
    nl2 = nls(k+1);
    x1  = (k-1);
    x2  = k;
    for j=1:nl1
        for m=1:nl2
            plot([x1 x2],[j-1-nl1/2 m-1-nl2/2] + bcent ,'k','linewidth',1);
        end
    end
end

for k = 1:numel(nls)
    nl1 = nls(k);
    y1  = (k-1);
    for j=1:nl1
        plot(y1,j-1 + bcent - nl1/2,'ko','markerfacecolor','w','color',call{k}(nl1+1-j,:),...
            'markersize',msize,'linewidth',4*msize/14);
    end
    for jj = 1:ells(k)
        plot(y1,-1*jj + bcent - nl1/2,'k.');
    end
end


axis tight;
axis off;
