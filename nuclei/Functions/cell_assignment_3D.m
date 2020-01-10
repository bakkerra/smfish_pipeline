function[cell_to_mRNA_map] = cell_assignment_3D(mRNA, vorvx)
% for each centroid, use the vertices from the voronoi diagram to determine
% which cell it belongs in. assign that centroid to the appropriate cell

fprintf('Assigning mRNA to cells\n');

% intialize structure to hold mRNA associated with each cell
% confusingly, use the matlab cell data structure
cell_to_mRNA_map = cell(length(vorvx), 1);

% initialize storage for each cell
% in the worst case, all of the centroids are in one cell
ZPlane_mod=mRNA.ZPlane*4.4;
centroids_3D=[mRNA.Centroid ZPlane_mod];
num_mRNA_found = zeros(length(vorvx), 1);

%make a 3d hull for each cell,find and list coordinates of mRNAs inside
%hull
for ii = 1:length(vorvx)
    bnd0 = vorvx{1,ii};
    in=inhull_mod(centroids_3D,bnd0);
    pnts=centroids_3D(in,:);
    num_mRNA_found(ii)=length(pnts);
    cell_to_mRNA_map{ii,1}=pnts; 
end

end