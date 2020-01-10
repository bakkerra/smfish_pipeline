function [mRNA] = driver(images_folder, thresh)
% pipeline to quantify the number of mRNA per cell in a volume of tissue

% inputs:
    % images_folder
        % string
        % relative path to folder contatining set of original 2D scans of 
        % fluorescent mRNA probes
    % [optional] thresh
        % integer
        % pixel value to use to threshold images into spots and background
% outputs
    % mRNA
        % struct
        % contains information about the n spots identified as mRNA
        % fields:
            % mRNA.Centroid
                % n by 2 double
                % each row [x y] is coordinates of centroids
            % mRNA.PixelIdxList
                % n by 1 cell
                % each row is array of pixels identified spot contains
            % mRNA.ZPlane
                % n by 1 double
                % each row is the z plane of that spots centroid
            % mRNA.IntensityAvg
                % n by 1 double
                % each row is the average pixel intensity of the pixels
                % within a circle of radius 4 around the centroid
            % mRNA.IntensityAvgCorrected
                % n by 1 double
                % mRNA.IntensityAvg with linear correction factor applied
                % to correct for bleaching
% example usage:
    % without optional threshold parameter:
        % [mRNA] = driver('Data/New_Probes/ATTO633_02/Images');
    % with optional threshold parameter:
        % [mRNA] = driver('Data/New_Probes/ATTO633_02/Images', 50);
        
% start timer
tic

% configure functions and libraries
config();

% fit line to the decrease in mean pixel intensity over z planes
[slope, images] = fit_bleaching(images_folder);

% if optional thresh parameter not provided, get from user
if nargin < 2
    % segment images at multiple thresholds, plot number of spots vs threshold
    sweep_thresholds(images, slope);
    thresh = input('Enter threshold to use to segment spots: ');
end

% segments images based on threshold, get centroids of spots
[candidates, num_candidates] = segmentation_mRNA(images, thresh);

% filter spots to determine which represent mRNA
[mRNA] = filter_mRNA(candidates, num_candidates, images, true);

% determine pixel intensities for each of the mRNA spots
[mRNA] = calculate_intensity(mRNA, images, slope);

% stop timer
toc
end
