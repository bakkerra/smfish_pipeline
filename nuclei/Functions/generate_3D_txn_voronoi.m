function [mRNA_nums]=generate_3D_txn_voronoi(vorvx,mRNA,image_size)
%same as generate_3D_txn_voronoi but generate a column that sums the
%nascent RNA in each cell.  2019/summed nascent rna not used for analysis
%for paper 


%segment nuclei
% [centroids_nuclei]=segmentation_nuclei(nuclei_h5);
% centroids_nuclei=centroid_pts;
%adjust z dimension to be properly equal to xy
% [centroids_nuclei(:,3)]=centroids_nuclei(:,3)*4;
%create 3D voronoi diagram 
% vorvx=bounded_3D_voronoi(centroids_nuclei,image_size);

fprintf('Assigning Txn Sites to Cells')
%assign txn mRNAs to voronoi "cells" in a Matlab cell
%data structure
mRNA_nums=struct;
[mRNA_nums.cells,mRNA_nums.totalnums]=cell_assignment_txn(mRNA,vorvx);

%calculate numbers of mRNAs per cell and place in a data structure
mRNA_nums.Totals=zeros(length(mRNA_nums.cells),1);
for i=1:length(mRNA_nums.cells)
    dime=size(mRNA_nums.cells{i,1});
    mRNA_nums.Totals(i,:)=dime(1);
    dime=[];
end

%calculate sum of txn mRNAs
for ii=1:length(mRNA_nums.totalnums)
    mRNA_nums.Sum(ii,:)=sum(mRNA_nums.totalnums{ii,1});
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