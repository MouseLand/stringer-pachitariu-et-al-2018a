% bins x in bins of tbin across idim dimension where idim = 1 or 2
function x = bin2d(x, tbin, idim)
if nargin<3
    idim = 1;
end   

if idim == 1
    x = reshape(mean(reshape(x(1:floor(size(x,1)/tbin)*tbin, :),...
        tbin,[],size(x,2)),1), [], size(x,2));
else
    x = squeeze(mean(reshape(x(:,1:floor(size(x,2)/tbin)*tbin),...
        size(x,1),tbin,[]),2));
end