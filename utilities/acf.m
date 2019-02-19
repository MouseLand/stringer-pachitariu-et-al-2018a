% autocorrelation of y (time x nd)
% L = number of timelags
function ac = acf(y, L)

mu = mean(y,1);
ac = NaN*zeros(L,size(y,2));
for k = 1:L
	ac(k,:) = gather_try(sum((y(k:end,:) - mu) .* (y(1:end-k+1,:) - mu),1));
end

yvar = sum((y - mu) .* (y - mu),1);
ac = gather_try(ac ./ yvar);
