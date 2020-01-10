function [centroids_nuclei,zplanes]=segmentation_driver(nuclei_tifs)
%generates a list of 3D nuclei centroids from segmented nuclei tifs
%imported from nucleiAIzer. Format for the tifs is 16bit imgs where each
%bit value corresponds to an 2D individual nuclei object.

%start timer
tic

%configure functions and libraries folder on path
config();

%connect overlapping 2d segmented nuclei
[bwnuclei,zplanes]=make_bw_images(nuclei_tifs);
[centroids_nuclei]=track_nuclei(bwnuclei,zplanes);

toc