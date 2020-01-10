function [mRNA_nums,txn_nums]=voronoi_driver(centroids_nuclei,zplanes,mRNA,txn_sites)
%Pipeline for generation of a 3D voronoi diagram from a set of nuclei
%centroids and 3d cell assignment of RNA and txn sites to cells. Returns
%'nums' data structures which are a list of cells

%convert zplanes to same scale as x-y pixels for voronoi
z_adj=ceil(zplanes*4.4);
image_size=[1024 1024 z_adj];
%total image size may need to be padded by adding 1-3 planes to prevent super small angles.
%no RNA spots will be in the 'padded' layers so doesn't affect cell assignment

%make 3D voronoi for regular spots
[mRNA_nums,vorvx]=generate_3D_voronoi(centroids_nuclei,mRNA,image_size);

%make 3D voronoi for transcription sites (recycles voronoi for speed)
[txn_nums]=generate_3D_txn_voronoi(vorvx,txn_sites,image_size);

%space coordinate
% [mRNA_nums]=space_coordinate(mRNA_nums,mRNA);
