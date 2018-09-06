# stringer-pachitariu-et-al-2018a

This code produces the figures from Stringer, Pachitariu et al, 2018a:

Carsen Stringer, Marius Pachitariu, Nicholas Steinmetz, Charu Bai Reddy, Matteo Carandini, Kenneth D. Harris
**Spontaneous behaviors drive multidimensional, brain-wide population activity**
https://www.biorxiv.org/content/early/2018/04/22/306019

It relies on data deposited on figshare at:

Carsen Stringer, Marius Pachitariu, Charu Bai Reddy, Matteo Carandini, Kenneth D. Harris
**Recordings of ten thousand neurons in visual cortex during spontaneous behaviors.** ([link](https://figshare.com/articles/Recordings_of_ten_thousand_neurons_in_visual_cortex_during_spontaneous_behaviors/6163622))

Here's the data [description](dataSharing.pdf).

The datasets to produce figures 1-4,6 are available. Shortly we will add new datasets for figure 5. 

The script 'mainFigs.m' calls all figure-producing scripts. Set useGPU=0 if you do not have a GPU.

While this is useful for reproducing the paper, it does not easily allow building on the main analysis we do, or running these analyses on different data. Therefore, we will add separate scripts to run analyses one at at time, with the data clearly loaded at the top of the script, so you can swap yours in. 

### How to load the data into python ###
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

<!---

Sorting by the top principal component.

### Analysis two ###

Reduced-rank regression from a set of behavioral variables to a set of neural activities.

### Analysis three ###

Low-dimensional nonlinear embedding into a 1D manifold. 

### Analysis four ###

Predicting a single neuron from its simultaneously-recorded peer using reduced-rank regression. 

### Analysis five ###

Comparison of population activity from two task conditions, for example spontaneous activity and stimulus presentation. 
 
--->
