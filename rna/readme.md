# RNA Segmentation

- To run an RNA segmentation, place a folder of your images in this folder
- The image folder should contain 16 bit black and white .tifs of a single channel of smFISH signal
- You must supply the title of this folder to the driver function
- You can optionally supply the threshold you want to segment at as a second argument.  If a threshold is not supplied, the function will segment at a range of thresholds and prompt your to pick one.  
- Generally, a good threshold is where there is an inflection point or plateau in the number of spots segmented at that channel
