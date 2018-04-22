function fx = acf_mat(x,L,varargin)
[NT] = size(x,1);
%x = single(x);
%x = (x - repmat(mean(x(:,:,:),2),1,size(x,2),1)) ./ ...
%repmat(std(x(:,:,:),[],2),1,size(x,2),1);
    %x = x/NT .^.5;
if isempty(varargin)
    x = zscore(x,1,1)/NT.^.5;
else
    x = x - repmat(mean(x,1),size(x,1),1);
end

fx = real(ifft(fft(x,[],1) .* fft(x(end:-1:1,:,:),[],1)));

fx = fx(1:round(NT/2),:,:);
fx = fx(1:L,:,:);

fx = gather(fx);
end