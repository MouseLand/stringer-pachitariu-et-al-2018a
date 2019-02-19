% where data is stored (that you download from figshare)
dataroot = '/media/carsen/DATA2/grive/10krecordings/spont_paper/';        
% where processed data and results are saved
matroot  = '/media/carsen/DATA2/grive/10krecordings/spontResults/';


%% this will produce the figures

close all;

%% fig1: 1st PC and SVCA
fig1new(matroot);

%% fig2: multi-dimensional behavior prediction in V1
fig2new(matroot);

%% fig3: ephys data
fig3new(matroot);

%% fig4: stim-spont comparison in V1
fig4new(matroot);

%% supplementary figures

% increasing neurons and bin size fig
suppIncNeuronsBins(matroot);

% drifting grating stim-spont
suppOriStimSpont(matroot);

% all the face embeddings
suppFacePreds(matroot);

% single neuron predictions
suppPeerNeurons(matroot);

% 1D behavioral variables correlations
supp1DBeh(matroot);

% PC autocorrelation functions
supppcs(matroot);

% distances between neurons vs correlations
suppDists(matroot);

% repeatability of correlations (1st half, 2nd half)
suppCorrStats(matroot);

% explainability by behavior vs position
suppVarPos(matroot);

% stim-spont variances
suppStimEmbed(matroot);

% ephys time lags
suppCCEphys(matroot);

% ephys manifold embedding separated by area
suppFastRasters(ephysroot, matroot)