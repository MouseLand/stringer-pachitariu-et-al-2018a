% cross-correlation of x and y (time x nd each)
% L = number of timelags
function cc = ccf(x, y, L)

mux = mean(x,1);
muy = mean(y,1);
cc = NaN*zeros(L*2+1,size(y,2));
for k = -L:L
	
	if k > 0
		cc(k+L+1,:) = gather_try(sum((y(k+1:end,:) - muy) .* (x(1:end-k,:) - mux),1));
	else
		cc(k+L+1,:) = gather_try(sum((y(1:end-abs(k),:) - muy) .* (x(abs(k)+1:end,:) - mux),1));
	end
end

yvar = sum((y - muy) .* (y - muy),1) .^ 0.5;
xvar = sum((x - mux) .* (x - mux),1) .^ 0.5;

cc = gather_try(cc ./ yvar ./ xvar );
