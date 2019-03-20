# stringer-pachitariu-et-al-2018a
[![DOI](https://zenodo.org/badge/130262484.svg)](https://zenodo.org/badge/latestdoi/130262484)


This code produces the figures from Stringer, Pachitariu et al, 2018a:

Carsen Stringer, Marius Pachitariu, Nicholas Steinmetz, Charu Bai Reddy, Matteo Carandini, Kenneth D. Harris
**Spontaneous behaviors drive multidimensional, brain-wide population activity**
https://www.biorxiv.org/content/10.1101/306019v2

It relies on data deposited on figshare at:

Carsen Stringer, Marius Pachitariu, Charu Bai Reddy, Matteo Carandini, Kenneth D. Harris
**Recordings of ten thousand neurons in visual cortex during spontaneous behaviors.** ([link](https://figshare.com/articles/Recordings_of_ten_thousand_neurons_in_visual_cortex_during_spontaneous_behaviors/6163622))

and

Nicholas Steinmetz, Marius Pachitariu, Carsen Stringer, Matteo Carandini, Kenneth D. Harris
**Eight probe neuropixels recordings during spontaneous behaviors.** ([link](https://figshare.com/articles/Eight-probe_Neuropixels_recordings_during_spontaneous_behaviors/7739750))

Here's the two-photon data [description](dataSharing.pdf). Please cite the paper and the figshare if you use the data.

The script 'processData.m' processes all the data. Set useGPU=0 if you do not have a GPU. The script 'mainFigs.m' calls all figure-producing scripts. 

While this is useful for reproducing the paper, it does not easily allow building on the main analysis we do, or running these analyses on different data. Therefore, we will add separate scripts to run analyses one at at time, with the data clearly loaded at the top of the script, so you can swap yours in. 

One such function that is easy to use is 'SVCA.m' which computes the cross-validated reliable variance of latent dimensions in the population activity.

### How to load the 2P data into python ###
```
import scipy.io as sio
import numpy as np
mt = sio.loadmat('spont_M150824_MP019_2016-04-05.mat')
spks = mt[‘Fsp’]                   # neurons by timepoints
med = mt[‘med’]                 # cell centers (X Y Z)

### behavioral measures
runSpeed = mt[‘beh’][‘runSpeed’]    # running speed
# get motion SVDs (time x components)
motionSVD=np.array(mt['beh'][0]['face'][0]['motionSVD'][0][0])  
# get motion masks (pixelY x pixelX x components)
motionMask=np.array(mt['beh'][0]['face'][0]['motionMask’][0][0])     
# pupil area and com
pupilArea =np.array(mt['beh'][0]['pupil'][0][‘area'][0][0])  
pupilCOM =np.array(mt['beh'][0]['pupil'][0][‘com'][0][0])  

# cell statistics
mt[‘stat’][0]     # first cell’s stats
mt[‘stat’][0][‘npix’]       # one example field, tells you how pixels make up the cell
```

### How to load the ephys data into python ###
```
from scipy import io

probeLoc = io.loadmat('/home/carsen/dm11/data/Spikes/eightprobes/probeLocations.mat')
probeBorders = io.loadmat('/home/carsen/dm11/data/Spikes/eightprobes/probeBorders.mat', squeeze_me=True)

mouse_names = ['Krebs','Waksman','Robbins']
# start of spontaneous activity in each mouse (in seconds)
tstart = [3811 3633 3323];

imouse = 0

spks = io.loadmat('/home/carsen/dm11/data/Spikes/eightprobes/spks/spks%s_Feb18.mat'%mouse_names[imouse], squeeze_me=True)
faces = io.loadmat('/home/carsen/dm11/data/Spikes/eightprobes/faces/%s_face_proc.mat'%mouse_names[imouse], squeeze_me=True)

# probe k
k = 0
# spike times (in seconds)
st = spks['spks'][k][0]
# clusters
clu = spks['spks'][k][1]
# cluster heights (in microns)
# (see siteCoords to convert to site location)
Wh = spks['spks'][k][2]

# processed faces
motSVD = faces['motionSVD']
video_timestamps = faces['times']

# where is the probe in the brain (consolidated labels)
# borders are in microns
# use Wh to determine which clusters are in which brain region
borders = probeBorders['probeBorders'][imouse]['borders'][k]
for j in range(len(borders)):
   b = borders[j]
   print('upper border %d, lower border %d, area %s'%(b[0],b[1],b[2]))
   wneurons = np.logical_and(Wh>=b[1], Wh<b[0])
   nn = wneurons.sum()
   print('%d neurons in %s'%(nn,b[-1]))
   
# where is the probe in the brain (in microns)
ccfCoords = probeLoc['probeLocations'][0][imouse]['probe'][k][0]['ccfCoords']
# name of area in Allen ontology by site on electrode
ccfNames = probeLoc['probeLocations'][0][imouse]['probe'][k][0]['ccfOntology']
# coordinates of each site on the electrode
siteCoords = probeLoc['probeLocations'][0][imouse]['probe'][k][0]['siteCoords']
```
