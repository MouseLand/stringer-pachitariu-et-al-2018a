% add the path of this folder
addpath(genpath('.'));

% you should change this to your local data paths
dataroot = '/media/carsen/DATA2/grive/10krecordings/spont_paper/';

% give a local folder for saving intermediate data (3GB max)
matroot = '/media/carsen/DATA2/mats/';

%% this will perform analyses and save output for figures

% figure 1
pcAnalysis(dataroot,matroot);

% figure 2 and 4

% 

%% this will produce the figures

fig1(matroot);

fig2(matroot);

fig3;

fig4;

%% supplementary figures