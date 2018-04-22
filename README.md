# stringer-pachitariu-et-al-2018a

This code produces the figures from Stringer, Pachitariu et al, 2018a:

Carsen Stringer, Marius Pachitariu, Nicholas Steinmetz, Charu Bai Reddy, Matteo Carandini, Kenneth D. Harris
**Spontaneous behaviors drive multidimensional, brain-wide population activity**
https://www.biorxiv.org/content/early/2018/04/22/306019

It relies on data deposited on figshare at:

Carsen Stringer, Marius Pachitariu, Charu Bai Reddy, Matteo Carandini, Kenneth D. Harris
**Recordings of ten thousand neurons in visual cortex during spontaneous behaviors. figshare. Fileset.**

Currently, only the main datasets are shared, which produce figures 1-4. Shortly we will add new datasets for figures 5 and 6. 

## Analyses ##

The script mainFigs calls all the analyses and the figure producing scripts. While this is useful for reproducing the paper, it does not easily allow building on the main types of analysis we do, or running them on different data. Therefore, we will add separate scripts to just run the following analyses, with the data clearly loaded at the top of the script, so you can swap yours in. 

### Analysis one ###

Sorting by the top principal component.

### Analysis two ###

Reduced-rank regression from a set of behavioral variables to a set of neural activities.

### Analysis three ###

Low-dimensional nonlinear embedding into a 1D manifold. 

### Analysis four ###

Predicting a single neuron from its simultaneously-recorded peer using reduced-rank regression. 

### Analysis five ###

Comparison of population activity from two task conditions, for example spontaneous activity and stimulus presentation. 
 