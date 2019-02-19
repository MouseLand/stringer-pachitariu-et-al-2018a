
function [v1,wpc,mot1,vc] = pc1Analysis(y,motSVD,tspont,tVid)

[NN NT] = size(y);
disp('>>>>>>>>>>>>> ');
disp([NN NT]);

%
ysub = mean(y,2);
y    = y - ysub;

% take SVDs of neural activity
[u s v] = svdecon(gpuArray(single(y)));
ncomps  = 128;
u       = gather(u);
s       = gather(s);
v       = gather(v);
u       = u(:, 1:ncomps);
v       = u' * y(:,:);
v1      = v(1,:);
wpc     = u(:,1);
% 
%% smooth motion energy PCs
z = motSVD;

%% correlation of first PC with arousal (motionSVD1)
td = -1;
tBeh  = td + tspont;
x = interp1(tVid, z, tBeh);
x = bin2d(x, 1.2*30, 1);

x1 = x(:,1);
x1 = sign(skewness(x1)) * x1;

%%

vc = corr(v(1,:)', x1);
wpc = sign(vc) * u(:,1);
v1 = sign(vc) * v1;
upc = sum(wpc>0)/NN;
disp([vc upc]);
mot1=x1;
%pc1 = v(1,:);

vc = corr(y',v1');

