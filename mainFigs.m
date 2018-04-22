addpath(genpath('.'));

% change to path of data
dataroot = '/media/carsen/DATA2/grive/10krecordings/spont_paper/';
% where the analyzed data saves
matroot = '/media/carsen/DATA2/mats/';

% do you have a GPU? if not set to 0
useGPU = 1;

%% this will perform analyses and save output for figures

%% figure 1
pcAnalysis(dataroot,matroot);

%% figure 2 and 4
predictNeuronsFromAllBeh(dataroot,matroot,useGPU);
quantifyBehavior(dataroot,matroot);

%% figure 3
smooth1Dclusters(dataroot,matroot,useGPU);
% peer prediction is pretty slow (different peers for each neuron)
% (I've put the mat file in the folder if you want to use it 'PCApred.mat')
peerExcludeNeighbors(dataroot,matroot,useGPU);

%% figure 4
% you can increase nseed in this script to average more but it will be
% slower (we used nseed = 10 in the paper)
faceTimelags(dataroot,matroot,useGPU);

%% this will produce the figures

fig1(matroot);

fig2(matroot);

fig3(matroot);
%%
fig4(matroot);

%% supplementary figures