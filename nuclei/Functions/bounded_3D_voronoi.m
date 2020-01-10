function vorvx = bounded_3D_voronoi(centroids_nuclei, image_size)
% build voronoi diagram of centroids taking into account the known
% boundaries of the image
 
fprintf('Building voronoi diagram\n')

% This code makes use of an open source project:
% Polytope bounded Voronoi diagram in 2D and 3D
% http://www.mathworks.com/matlabcentral/fileexchange/50772-polytope-bounded-voronoi-diagram-in-2d-and-3d

% Copyright (c) 2015, Hyongju Park
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% boundaries from size of image
bnd0 = [1 1 1; image_size(2) 1 1; 1 image_size(1) 1; image_size(2) image_size(1) 1; 1 1 image_size(3); image_size(2) 1 image_size(3); 1 image_size(1) image_size(3); image_size(2) image_size(1) image_size(3)];
K = convhull(bnd0);
bnd_pnts = bnd0(K,:);

in = inhull(centroids_nuclei,bnd0); 
% inhull.m is written by John D'Errico that efficiently check if points are
% inside a convex hull in n dimensions
% We use the function to choose points that are inside the defined boundary
u1 = 0;
for i = 1:size(centroids_nuclei,1)
    if in(i) ==1
        u1 = u1 + 1;
        pos(u1,:) = centroids_nuclei(i,:);
    end
end

% =========================================================================
% INPUTS:
% pos       points that are in the boundary      n x d matrix (n: number of points d: dimension) 
% bnd_pnts  points that defines the boundary     m x d matrix (m: number of vertices for the convex polytope
% boundary d: dimension)
% -------------------------------------------------------------------------
% OUTPUTS:
% vornb     Voronoi neighbors for each generator point:     n x 1 cells
% vorvx     Voronoi vertices for each generator point:      n x 1 cells
% =========================================================================

[~,vorvx] = polybnd_voronoi(pos,bnd_pnts);

% figure;
% for i = 1:size(pos,1)
%     plot(vorvx{i}(:,1),vorvx{i}(:,2),'-r')
%     hold on;
% end
% plot(bnd_pnts(:,1),bnd_pnts(:,2),'-');
% hold on;
% scatter(pos(:,1),pos(:,2),'Marker','o','MarkerFaceColor',[0 .75 .75],'MarkerEdgeColor','k');
% axis equal;
% axis([0 image_size(2) 0 image_size(1)]);
% set(gca,'YDir','reverse'); % flips y coordinates so plot matches image file

for i = 1:size(vorvx,2)
    col(i,:)= rand(1,3);
end

figure;
for i = 1:size(pos,1)
    K = convhulln(vorvx{i},{'QJ'});
    trisurf(K,vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'FaceColor',col(i,:),'FaceAlpha',0.5,'EdgeAlpha',1)
    hold on;
end
scatter3(pos(:,1),pos(:,2),pos(:,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
axis equal;
axis([0 image_size(2) 0 image_size(1) 0 image_size(3)]);
set(gca,'YDir','reverse'); % flips y coordinates so plot matches image file
xlabel('X');ylabel('Y');zlabel('Z');