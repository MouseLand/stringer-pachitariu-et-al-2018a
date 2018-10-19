clear all

% you should change this to your local data paths
dataroot = '/media/carsen/DATA2/grive/10krecordings/spontData';
%dataroot = 'D:/grive/10krecordings/spontData';

% give a local folder for saving intermediate data (3GB max)
matroot = '/media/carsen/DATA2/grive/10krecordings/spontResults';
%matroot = 'D:/grive/10krecordings/spontResults';

mkdir(matroot)

% do you have a GPU? if not set to 0
useGPU = 1;

% should be in github folder
addpath(genpath('.'));

% also download rastermap
% https://github.com/MouseLand/rastermap/
addpath('/media/carsen/DATA2/Github/rastermap/matlab/');
%addpath('C:\Users\carse\github\rastermap\matlab');

dex = 2; % second dataset as example dataset

%% this will perform analyses and save output for figures

%% run PC analysis
pcAnalysis(dataroot,matroot, useGPU, dex);

%% run rastermap and determine distance-dependence of clusters
smooth1Dclusters(dataroot,matroot,useGPU, dex);

%% predict PCs from each other
peerPC_CCA(dataroot, matroot, useGPU, dex);

%% peer prediction for neurons is pretty slow (different peers for each neuron)
% (I've put the mat file in the folder if you want to use it 'PCApred.mat')
% (this is only in a supplementary fig)
peerExcludeNeighbors(dataroot,matroot,useGPU);

%% run behavioral analyses
predictNeuronsFromAllBeh(dataroot,matroot,useGPU, dex);
quantifyBehavior(dataroot,matroot, dex);

%% time delay analysis (panel 4K)
% you can increase nseed in this script to average more but it will be
% slower (we used nseed = 10 in the paper)
faceTimelags(dataroot,matroot,useGPU);

%% run analysis for figure 6
sharedVariance(dataroot, matroot, useGPU);
faceStatistics(dataroot, matroot);
stimfaceVariance(dataroot,matroot,useGPU);

