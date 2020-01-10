function [slope, images] = fit_bleaching(images_folder)
% fit line to the decrease in mean pixel intensity over z planes

fprintf('Fitting bleaching curve\n');

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
mean_intens = zeros(1,size(images,3));
for kk = 1:size(images,3)
    mean_intens(kk) = mean2(images(:,:,kk));
end
figure;
hold on;
x = 1:length(mean_intens);
y = mean_intens;
polynomial_coeffs = polyfit(x,y,1);
fittedy = polyval(polynomial_coeffs,x);
r2 = rsquare(y,fittedy);
slope = polynomial_coeffs(1);
scatter(x, y, 'filled');
plot(x, fittedy);
title('Bleaching Curve');
ylabel('Mean Intensity')
xlabel('Z Plane');
ylim([0 inf]);
text(2,2, sprintf('Slope: %.3f\nr^2: %.3f', slope, r2));

cd(old_path);