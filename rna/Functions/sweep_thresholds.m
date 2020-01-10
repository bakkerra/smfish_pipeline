function sweep_thresholds(images, slope)
% segment images at multiple thresholds, plot number of spots vs threshold

fprintf('Sweeping Thresholds \nThis could take a little while...\ngrab a cup of coffee and remember the threshold you chose at the end\nso if you run this dataset again you can supply the threshold as an argument and not wait\n');
% apply linear correction factor for bleaching
corrected_images = uint16(zeros(size(images)));
for ii = 1:size(images,3)
    corrected_images(:,:,ii) = uint16(images(:,:,ii) - slope * ii);
end

threshs = 15:150;%1:max(max(max(corrected_images)));

% calculate number of spots identified at each threshold
num_identifieds_orig = zeros(size(threshs));
num_identifieds_corrected = zeros(size(threshs));
for kk = 1:length(threshs)
    fprintf('Thresh: %d ', kk+25-1);
    thresh = threshs(kk);
    identified_centroids_orig = identify_centroids_thresh(images, thresh);
    identified_centroids_corrected = identify_centroids_thresh(corrected_images, thresh);
    
    num_identifieds_orig(kk) = length(identified_centroids_orig);
    num_identifieds_corrected(kk) = length(identified_centroids_corrected);
end

% plot should plateau at ideal threshold 
figure;
hold on;
plot(threshs,num_identifieds_orig);
plot(threshs,num_identifieds_corrected);
title('Number of Spots Identified vs Threshold');
xlabel('Threshold Parameter (Pixel Value)');
ylabel('Number of Spots Identified');
legend('Original', 'Bleaching Corrected')
hold off;
end

function [candidates, num_candidates, bw_thresholded_mRNA, image_size] = threshold_raw_images(images, thresh)
    image_size = size(images);
    bw_thresholded_mRNA = false(image_size);
    for ii=1:image_size(3)
        bw_thresholded_mRNA(:,:,ii) = im2bw(images(:,:,ii), thresh);
        % find connected components, identify their centroids and list of pixels
        CC = bwconncomp(bw_thresholded_mRNA(:,:,ii),8); %connectivity 8 in 2D
        stats = regionprops(CC, 'Centroid', 'PixelIdxList');
        data_mRNA(ii).stats = stats;
    end

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
end

function identified_centroids = identify_centroids_thresh(images, thresh_pix)
    thresh = thresh_pix/2^16; % images are 16 bit
    [candidates, num_candidates] = threshold_raw_images(images, thresh);
    [mRNA] = filter_mRNA(candidates, num_candidates, images, false);
    identified_centroids = [mRNA.Centroid(:,2), mRNA.Centroid(:,1),  mRNA.ZPlane];
end