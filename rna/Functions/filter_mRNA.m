function[mRNA] = filter_mRNA(candidates, num_candidates, images, DISPLAY)
% filter spots to determine which represent mRNA

fprintf('Filtering mRNA spots\n');

MIN_SIZE = 5; % pixels

% go through each plane. for each spot in that plane,
% keep only ones with at least 1 neighbor in an adjacent plane
% and that are have the maximum size of all their neighbors
good_centroids = zeros(num_candidates, 2);
good_pixels = cell(num_candidates,1);
good_zplanes = zeros(num_candidates,1);
num_good = 0;
for ii = 1:length(candidates)
    for jj = 1:size(candidates(ii).Centroid,1)
        bot_size = 0; % default values so still works on edge planes
        top_size = 0;
        bot_found = false;
        top_found = false;
        
        if ii == 1
            bot_found = true;
        else % if not bottom plane
            % check lower plane
            bot_size = adjacent_size(candidates(ii).Centroid(jj,:), ii - 1);
            if bot_size ~= 0
                bot_found = true;
            end
        end
        
        if ii == length(candidates)
            top_found = true;
        else % if not top plane
            % check upper plane
            top_size = adjacent_size(candidates(ii).Centroid(jj,:), ii + 1);
            if top_size ~= 0
                top_found = true;
            end
        end
        
        if bot_found || top_found
            % if has at least one neighbor and has max spot size of neighbors
            curr_size = length(candidates(ii).PixelIdxList{jj, 1});
            if (curr_size > bot_size) && (curr_size > top_size) && (curr_size > MIN_SIZE)
                % doesn't count at all in case of equal sizes in adjacent spots
                num_good = num_good + 1;
                good_centroids(num_good,:) = candidates(ii).Centroid(jj,:);
                good_pixels{num_good,1} = candidates(ii).PixelIdxList{jj,1};
                good_zplanes(num_good) = ii;
            end
        end
    end
end

% get rid of empty pre allocated spaces and put in final data structure
% store everything in the data structure for mRNAs.
mRNA = struct;
mRNA.Centroid = good_centroids(1:num_good,:);
mRNA.PixelIdxList = good_pixels(1:num_good,1);
mRNA.ZPlane = good_zplanes(1:num_good);

    % Note: nested function, has access to the "global" variable candidates
    function max_size = adjacent_size(centroid, adj_ind)
        % determines if a centroid has a neighbor in a given adjacent plane
        
        % inputs: centroid is [x y]
        %         adj_ind is the index of the plane searching for neighbors in
        % output: returns size of spot if finds neighbor, else returns zero
        %         if finds more than one neighbor, returns size of largest
        
        % centroids must less than MAX_DIST pixels appart to be considered neighbors
        MAX_DIST = 4; 
        
        curr_x = centroid(1);
        curr_y = centroid(2);
        
        % mask centroids within square to get rough estimate of dist
        % don't have to check as many points
        x_mask = (candidates(adj_ind).Centroid(:,1) <= curr_x + MAX_DIST) ...
            & (candidates(adj_ind).Centroid(:,1) >= curr_x - MAX_DIST);
        y_mask = (candidates(adj_ind).Centroid(:,2) <= curr_y + MAX_DIST) ...
            & (candidates(adj_ind).Centroid(:,2) >= curr_y - MAX_DIST);
        mask = x_mask & y_mask;
        pretty_close = candidates(adj_ind).Centroid(mask,:);
        pretty_close_pix_lists = candidates(adj_ind).PixelIdxList(mask);
        
        % check if any centroids are actually within allowable distance
        % if so, return size of largest spot
        max_size = 0;
        for nn = 1:size(pretty_close,1)
            if sqrt( (curr_x - pretty_close(nn,1))^2 + (curr_y - pretty_close(nn,2))^2 ) <= MAX_DIST
                spot_size = length(pretty_close_pix_lists{nn, 1});
                if spot_size > max_size
                    max_size = spot_size;
                end
            end
        end
    end
if DISPLAY
    size_hist(mRNA);
    %reconstruct_image(mRNA,image_size);
    dot_images(mRNA, images);
end
end

function size_hist(mRNA)
% generates histogram of frequency of sizes of segmented objects
    n = length(mRNA.PixelIdxList);
    sizes = zeros(n,1);
    for ii = 1:n
        sizes(ii) = size(mRNA.PixelIdxList{ii,1},1);
    end
    figure, histogram(sizes); %unfiltered
    title('Sizes of Identified Spots (After Filtering)');
    xlabel('Size of Spot (pixels)')
    ylabel('Number of Spots');
    
% uncomment to view filtered histograms    
%     figure, histogram(sizes(sizes < 500));
%     title('Frequency of spot sizes less than 500 pixels');
%     xlabel('Size of spot (pixels)')
%     ylabel('Number of spots');
%     
%     figure, histogram(sizes(sizes < 150));
%     title('Frequency of spot sizes less than 150 pixels');
%     xlabel('Size of spot (pixels)')
%     ylabel('Number of spots');
end

function dot_images(mRNA, images)
    centroids = mRNA.Centroid;
    images_8bit = uint8(images); % matlab can only display 8 bit images
    dotted_image = uint8(zeros(size(images_8bit,1), size(images_8bit,2), 3, size(images_8bit,3)));
    marker = false(size(images_8bit));
    for k = 1:length(centroids)
        xy = round(centroids(k,:));
        marker(xy(2), xy(1), mRNA.ZPlane(k)) = 1;
    end
    
    % make dots larger
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

