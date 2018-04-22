addpath(genpath('.'));

% change to path of data
dataroot = '/media/carsen/DATA2/grive/10krecordings/spont_paper/';
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