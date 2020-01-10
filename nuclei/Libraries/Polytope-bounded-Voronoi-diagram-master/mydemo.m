clear all;close all;clc

%% generate random samples
n = 5;        % number of points
m = 4;         % number of boundary point-candidates
d = 2;          % dimension of the space
pos0 = [50,50; 40,900; 250,500; 480,60];       
bnd0 = [1 1; 512 1; 1 1024; 512 1024];       % generate boundary point-candidates
K = convhull(bnd0);
bnd_pnts = bnd0(K,:);   % take boundary points from vertices of convex polytope formed with the boundary point-candidates
%% take points that are in the boundary convex polytope
in = inhull(pos0,bnd0); 
% inhull.m is written by John D'Errico that efficiently check if points are
% inside a convex hull in n dimensions
% We use the function to choose points that are inside the defined boundary
u1 = 0;
for i = 1:size(pos0,1)
    if in(i) ==1
        u1 = u1 + 1;
        pos(u1,:) = pos0(i,:);
    end
end
%% 
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

[vornb,vorvx] = polybnd_voronoi(pos,bnd_pnts);

%% PLOT




figure;
for i = 1:size(pos,1)
plot(vorvx{i}(:,1),vorvx{i}(:,2),'-r')
hold on;
end
plot(bnd_pnts(:,1),bnd_pnts(:,2),'-');
hold on;
scatter(pos(:,1),pos(:,2),'Marker','o','MarkerFaceColor',[0 .75 .75],'MarkerEdgeColor','k');
axis equal;
axis([0 512 0 1024]);
set(gca,'YDir','reverse'); % flips y coordinates so plot matches image file      

