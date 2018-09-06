% add the code repository to the path
addpath(genpath('D:\Github\stringer-pachitariu-et-al-2018a\'))

% this loads my compiled data file, that contains "Ff" = neurons by time
load('H:\spont_MP019.mat')

useGPU = 0; % set this to 1 if you have an Nvidia GPU set up in Matlab

nC = 15; % number of clusters along 1D embedding: it can get confused if this is too large. 

%% 
S = zscore(Ff, 1, 2)/size(Ff,2).^.5;
[NN, NT] = size(S);

if useGPU
    S = gpuArray(single(S));
end

% the top PC is used to initialize the ordering
[U, ~,~] = svdecon(S); % svdecon is contained in the repository
U        = U(:,1); 
[~, isort] = sort(U(:,1)); 

% run the embedding
[iclustup, isort] = embed1D(S, nC, isort, useGPU); % 

%%  this cell plots the cells sorted by the ordering and smoothed over cells
Sm = S;
Sm = Sm(isort, :);

if useGPU
    Sm = gpuArray(Sm);
end

sig_cells = 10; % replace this with whatever value makes the plot look good! 
Sm = my_conv2(Sm, sig_cells, 1); 

imagesc(Sm, [-.2 .7]/150) % play with the scaling of the image too
cmap = colormap('gray');
colormap(flipud(cmap))

%%