function bin_mRNA(mRNA, image_size)
% bin the mRNA into stripes and grids to compare their varience
% (as opposed to assigning to cells)

fprintf('Binning mRNA THIS STILL NEEDS TO BE DEVELOPED FULLY\n');

NUM_X_BINS = 32;
NUM_Y_BINS = 32;

x_bin_width = round(image_size(2) / NUM_X_BINS);
y_bin_width = round(image_size(1) / NUM_Y_BINS);

num_spots_x = zeros(NUM_X_BINS, 1);
num_spots_y = zeros(NUM_Y_BINS, 1);
num_spots_grid = zeros(NUM_Y_BINS, NUM_X_BINS);
for ii = 1:size(mRNA.Centroid, 1)
    x_coord = mRNA.Centroid(ii, 1);
    x_bin = ceil(x_coord/x_bin_width);
    if x_bin <= NUM_X_BINS
        num_spots_x(x_bin) = num_spots_x(x_bin) + 1;
    else
        printf('fix this case\n');
    end
    
    y_coord = mRNA.Centroid(ii,2);
    y_bin = ceil(y_coord/y_bin_width);
    if y_bin <= NUM_Y_BINS
        num_spots_y(y_bin) = num_spots_y(y_bin) + 1;
    else
        printf('fix this case\n');
    end
    
    if (y_bin <= NUM_Y_BINS) && (x_bin <= NUM_X_BINS)
        num_spots_grid(y_bin, x_bin) = num_spots_grid(y_bin, x_bin) + 1;
    else
        printf('fix this case\n')
    end
    
end
figure, plot(1:NUM_X_BINS, num_spots_x);
title('Number of Spots in DV Bins of Width 16 px');
xlabel('Bin in DV (Width 16 px)');
ylabel('Number of spots');
ylim([0 inf]);
xlim([1 NUM_X_BINS]);

figure, plot(1:NUM_Y_BINS, num_spots_y);
title('Number of Spots in AP Bins of Width 16 px');
xlabel('Bin in AP (Width 16 px)');
ylabel('Number of spots');
ylim([0 inf]);
xlim([1 NUM_Y_BINS]);

figure, surf(num_spots_grid);
title('Number of Spots in 16 px x 16px Bins')
xlabel('Bin in DV (Width 16 px)');
ylabel('Bin in AP (Width 16 px)');
zlabel('Number of spots');
c = colorbar;
c.Label.String = 'Number of Spots per Bin';
axis equal
axis([1 NUM_X_BINS 1 NUM_Y_BINS]);

