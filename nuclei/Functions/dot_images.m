function dot_images(mRNA, images_folder)
    
old_path = cd(images_folder);

image_dir = dir;
num_images = 0;
for jj = 1:length(image_dir)
    if image_dir(jj).name(1) ~= '.' % ignore hidden files
        num_images = num_images + 1;
        flipped_image = imread(image_dir(jj).name);
        % transpose to orient original images the same way as the h5
        images(:,:,num_images) = flipped_image';
    end
end


centroids = mRNA.Centroid;
images_8bit = uint8(images); % matlab can only display 8 bit images
dotted_image = uint8(zeros(size(images_8bit,1), size(images_8bit,2), 3, size(images_8bit,3)));
marker = false(size(images_8bit));
for k = 1:length(centroids)
    xy = round(centroids(k,:));
    marker(xy(2), xy(1), mRNA.ZPlane(k)) = 1;
end
% make dots larger
%disk2 = strel('disk',2);
%marker = imdilate(marker, disk2);
    
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
cd(old_path);
    