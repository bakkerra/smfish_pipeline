function [bwnuclei,zplanes]=make_bw_images(nuclei_tifs)
%uses segmented images from nucleAIzer to make bw image

% import labelled nuclei images preprocessed from nucleAIzer
old_path = cd(nuclei_tifs);
image_dir = dir;
num_images = 0;
for jj = 1:length(image_dir)
    if image_dir(jj).name(1) ~= '.' % ignore hidden files
        num_images = num_images + 1;
        flipped_image = imread(image_dir(jj).name);
        % transpose to orient original images the same way as the h5
        %4/2018 this is outdated now, but FISH data is already analyzed
        %transposed
        images(:,:,num_images) = flipped_image';
    end    
end
zplanes=num_images;
cd(old_path);

%convert greyscale to bw 
se=strel('disk',1);
bwnuclei=false([1024,1024,zplanes]);
for i=1:zplanes
    current_img=images(:,:,i);
    objects=unique(current_img);
    current_final=false([1024 1024]);
    for j=1:max(objects)
        current_obj=current_img==j;
        shrinked_obj=imerode(current_obj,se);
        current_final(shrinked_obj)= true;
    end
    bwnuclei(:,:,i)=current_final;
end
%display results
implay(bwnuclei)
    