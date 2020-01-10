function[mRNA] = calculate_intensity(mRNA,images,slope)
% determine pixel intensity of each of the mRNA spots by sampling original
% images in locations determined to represent mRNA

fprintf('Determining mRNA pixel intensities\nProgress:\n');

% mask to visualize which pixels are selected
mask = false(size(images));

% build a circular structuring element (se) to be used around each centroid
% to sample pixels
se = strel('Disk', 4);

% sample original images to determine pixel intensity values for each spot.
intensity_avg = zeros(size(mRNA.Centroid,1), 1);
intensity_avg_corrected = zeros(size(mRNA.Centroid,1), 1);
curr_plane = 1;
curr_image = images(:,:,curr_plane);
for ii = 1:size(mRNA.Centroid, 1)
    %print out iteration number occasionally to keep track of where you are
    %in the function
    if (mod(ii,100) == 0)
        fprintf('%i/%i\n', ii, size(mRNA.Centroid,1)) 
    end
    
    % since we're using circle (2D) as structuring element,
    % only use 2D images (smaller) to minimize number of disk references
    % (slow)
    if mRNA.ZPlane(ii) > curr_plane
        curr_plane = mRNA.ZPlane(ii);
        curr_image = images(:,:,curr_plane);
    end
    
    % centroid is not always at a round number because of the averaging of
    % the connected components
    centroid = floor(mRNA.Centroid(ii,:));
    for j = 1:length(centroid)
        if centroid(j) == 0
            centroid(j) = 1;
        end
    end
    
    % build a binary image with just a point at the given centroid
    find_pixels = false(size(curr_image));
    find_pixels(centroid(2), centroid(1)) = true;
    
    % use the structuring element to dilate the centroid pixel into desired
    % shape. add to mask to visualize selected pixels
    find_pixels_dil = imdilate(find_pixels,se);
    mask(:,:,curr_plane) = mask(:,:,curr_plane) | find_pixels_dil; % slow
    
    % select pixels
    pixel_intensities = curr_image(find_pixels_dil);

    %calculate average and total intensity
    intensity_avg(ii) = mean(pixel_intensities);
    %intensity_tot(ii) = sum(pixel_intensities); really measuring the same
    % thing since the number of pixels sampled is always the same
    % i'll go with avg to make it easier to correct using the bleaching
    % slope
    
    % correct for bleaching
    intensity_avg_corrected(ii) = intensity_avg(ii) - slope * curr_plane;
end

implay(mask);

% generate histograms of brightness of spots
figure, histogram(intensity_avg);
title('Average Intensity of Spots');
xlabel('Average Pixel Intensity');
ylabel('Number of Spots');
figure, histogram(intensity_avg_corrected);
title('Average Intensity of Spots (Bleaching Corrected)');
xlabel('Corrected Average Pixel Intensity');
ylabel('Number of Spots');


%save it in the data structure for mRNA properties
mRNA.IntensityAvg = intensity_avg;
mRNA.IntensityAvgCorrected = intensity_avg_corrected;
end
