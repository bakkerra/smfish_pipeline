function [mRNA_nums,vorvx]=generate_3D_voronoi(centroids_nuclei,mRNA,image_size)

% segment nuclei (uncomment if using outside driver)
% [centroids_nuclei]=segmentation_nuclei(nuclei_h5);

fprintf('Making Voronoi Diagram\n')

%initialize data structure
mRNA_nums=struct;
mRNA_nums.Centroid=centroids_nuclei;

%adjust z dimension to be properly equal to xy
[centroids_nuclei(:,3)]=centroids_nuclei(:,3)*4.4;

%create 3D voronoi diagram 
vorvx=bounded_3D_voronoi(centroids_nuclei,image_size);

fprintf('Assigning mRNAs to Cells\n')

%assign mRNAs to voronoi "cells" in a Matlab cell
%data structure
[mRNA_nums.cells]=cell_assignment_3D(mRNA,vorvx);

%calculate numbers of mRNAs per cell and place in a data structure
mRNA_nums.Totals=zeros(length(mRNA_nums.cells),1);
for i=1:length(mRNA_nums.cells)
    dime=size(mRNA_nums.cells{i,1});
    mRNA_nums.Totals(i,:)=dime(1);
    dime=[];
end

%display figures with voronoi cells shaded with numbers
% axis_size = [0, image_size(2), 0, image_size(1)];
% figure;
% for jj = 1:size(vorvx,2)
%     patch(vorvx{1,jj}(:,1),vorvx{1,jj}(:,2), mRNA_nums.Totals(jj));
%     set(gca,'YDir','reverse');
%     axis equal
%     axis(axis_size)
%     c = colorbar;
%     c.Label.String = 'mRNA per Cell';
% end

figure;
for i = 1:size(vorvx,2)
    K = convhulln(vorvx{i},{'QJ'});
    trisurf(K,vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),mRNA_nums.Totals(i),'FaceAlpha',0.5,'EdgeAlpha',1)
    hold on;
end
axis equal;
axis([0 image_size(2) 0 image_size(1) 0 image_size(3)]);
c=colorbar;
set(gca,'YDir','reverse'); % flips y coordinates so plot matches image file
xlabel('X');ylabel('Y');zlabel('Z');