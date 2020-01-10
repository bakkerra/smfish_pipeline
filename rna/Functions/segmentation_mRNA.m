function [candidates, num_candidates] = segmentation_mRNA(images, thresh)
% segments images based on threshold, get centroids of spots

fprintf('Segmenting mRNA spots\n');

image_size = size(images);
% segment each z slice individually
num_objects_per_plane = zeros(image_size(3),1);
bw_thresholded_mRNA = false(image_size);

for ii=1:image_size(3)
    %thresh = graythresh(prob_image_mRNA(:,:,ii));
    % thresh is manually selected now, and has to be converted to fraction
    frac_thresh = thresh/2^16; %16 bit images
    bw_thresholded_mRNA(:,:,ii) = im2bw(images(:,:,ii), frac_thresh);
    
    % find connected components, identify their centroids and list of pixels
    CC = bwconncomp(bw_thresholded_mRNA(:,:,ii),8); %connectivity 8 in 2D
    stats = regionprops(CC, 'Centroid', 'PixelIdxList');
    data_mRNA(ii).stats = stats;
    num_objects_per_plane(ii) = CC.NumObjects;
end

% graph number of objects segmented per z plane
figure, scatter(1:length(num_objects_per_plane), num_objects_per_plane, 'filled');
title('Number of Objects Segmented vs. Z Plane');
xlabel('Z Plane');
ylabel('Number of Objects Segmented');
ylim([0 inf]); % force y axis to start at 0

% display segmented image
reconstruct_image(data_mRNA, image_size);

% convert centroid and pixel stats into easier to use datatype, candidates.
% will use in threshold to select which candidates represent mRNA
candidates = struct;
num_candidates = 0;
for ii = 1:length(data_mRNA) % for each z plane in image
    num_centroids = length(data_mRNA(ii).stats);
    centroid = zeros(num_centroids, 2);
    pixels = cell(num_centroids, 1);
    for jj = 1:num_centroids % for each centroid in that plane
        centroid(jj,:) = [data_mRNA(ii).stats(jj).Centroid];
        pixels{jj,1} = data_mRNA(ii).stats(jj).PixelIdxList;
    end
    candidates(ii).Centroid = centroid;
    candidates(ii).PixelIdxList = pixels;
    
    % keep track of total number of candidates found across all planes
    num_candidates = num_candidates + num_centroids;
end
dot_images(candidates, images);
size_hist(candidates);
end

function reconstruct_image(data_mRNA, image_size)
    image = zeros(image_size);
    for ii = 1:length(data_mRNA)
        plane_correction = (ii - 1) * image_size(1) * image_size(2);
        for jj = 1:length(data_mRNA(ii).stats)
            image(data_mRNA(ii).stats(jj).PixelIdxList + plane_correction) = 1;
        end
    end
    implay(image);
%    save video
%    v = VideoWriter('newfile.avi');
%    open(v);
%    for kk = 1:size(image,3)
%        writeVideo(v,image(:,:,kk));
%    end
%    close(v);
end  

function dot_images(candidates, images)
    images_8bit = uint8(images); % matlab can only display 8 bit images
    dotted_image = uint8(zeros(size(images_8bit,1), size(images_8bit,2), 3, size(images_8bit,3)));
    marker = false(size(images_8bit));
    for k = 1:length(candidates)
        for j = 1:size(candidates(k).Centroid, 1)
            xy = round(candidates(k).Centroid(j,:));
            marker(xy(2), xy(1), k) = 1;
        end
    end
    
%     make dots larger
    disk2 = strel('disk',2);
    marker = imdilate(marker, disk2);
    
    rgb = [2^8, 0, 0];
    red = images_8bit;
    red(marker) = rgb(1);
    green = images_8bit;
    green(marker) = rgb(2);
    blue = images_8bit;
    blue(marker) = rgb(3);
    for i = 1:size(images_8bit,3)
        dotted_image(:,:,1,i) = red(:,:,i);
        dotted_image(:,:,2,i) = green(:,:,i);
        dotted_image(:,:,3,i) = blue(:,:,i);
    end
    implay(dotted_image)
    
end

function size_hist(mRNA)
% generates histogram of frequency of sizes of segmented objects
    sizes = [];
    for jj = 1:size(mRNA,2)
        inside_sizes = zeros(length(mRNA(jj).PixelIdxList), 1);
        for ii = 1:length(mRNA(jj).PixelIdxList)
            inside_sizes(ii) = size(mRNA(jj).PixelIdxList{ii,1},1);
        end
        sizes = [sizes; inside_sizes];
    end
        figure, histogram(sizes); %unfiltered
        title('Sizes of Candidate Objects Segmented (Before Filtering)');
        xlabel('Size of Object (pixels)')
        ylabel('Number of Objects');
    
end
    