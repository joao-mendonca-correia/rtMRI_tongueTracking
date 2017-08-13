# rtMRI_tongueTracking
Processes magnetic resonance (MR) images obtained using real-time MRI protocols during speech production, and tracks movements of the tongue based on small initial user input.

The aim of this method is to process a sequence of rtMRI images, frame by frame, into a set of isolated large clusters using image binarization. The cluster with the higher level of similarity to the tongue cluster estimated in the previous frame is selected as the current tongue cluster. The procedure continuous recursively. The first tongue cluster is obtained from user selection using mouse clicks.

To avoid the inclusion of other vocal-tract structures into the estimation of the tongue cluster, a dynamic localization of the lower and upper lips is done. Using a given lip thickness in pixels, the algorithm restricts the possibility of tongue clusters to invade lip areas.

The user is also allowed to indicate an area for which tongue cluster estimation is not allowed to enter, which may be useful to avoid expansion onto palate and velum regions.

The user must initialize the algorithm in 3 steps:
- Step 1: select a squared region of interest for which tongue segmentation will be performed. 
- Step 2: segment the initial tongue shape from the first frame using mouse clicks.
- Step 3: indicate a region for which recursive tongue segmentation should not enter using mouse clicks.

After this initialization, the algorithm computes tongue segmentation for every frame in the input path, as well as, lower-lips (LL) and upper-lips (UL).

Output files:

- outputMasks.mat:
    - tongue mask per frame
    - lower-lip mask per frame
    - upper-lip mask per frame
    - estimation of lip apperture
    
- outputMovie.mat:
  - Movie
 
The movie can be displayed using the command: implay(Movie);
