function disp_results(vorvx, cell_to_mRNA_map, mRNA, image_size)
% generate heatmaps displaying results

fprintf('Visualizing results\n');

% prepare data for display
num_mRNA = zeros(length(cell_to_mRNA_map));
intens_mRNA = zeros(length(cell_to_mRNA_map));
area_cell = zeros(length(cell_to_mRNA_map));
for ii = 1:length(cell_to_mRNA_map); % for each cell
    % determine number of mRNA in the cell
    num_mRNA(ii) = size(cell_to_mRNA_map{ii,1},2);
    
    % sum intensity of mRNA in the cell
    intens_sum = 0;
    for kk = 1:size(cell_to_mRNA_map{ii,1},2)
        intens_sum = intens_sum + mRNA.IntensityTot(1,cell_to_mRNA_map{ii,1}(1,kk));
    end
    intens_mRNA(ii) = intens_sum;
    
    % calculate area of the cell
    area_cell(ii) = polyarea(vorvx{1,ii}(:,1),vorvx{1,ii}(:,2));
end

% create figures
axis_size = [0, image_size(1), 0, image_size(2)];
figure;
for jj = 1:size(vorvx,2)
    patch(vorvx{1,jj}(:,1),vorvx{1,jj}(:,2), num_mRNA(jj));
    set(gca,'YDir','reverse');
    axis equal
    axis(axis_size)
    c = colorbar;
    c.Label.String = 'Number of mRNA per Cell';
end

figure;
for jj = 1:size(vorvx,2)
    patch(vorvx{1,jj}(:,1),vorvx{1,jj}(:,2), intens_mRNA(jj));
    set(gca,'YDir','reverse');
    axis equal
    axis(axis_size)
    c = colorbar;
    c.Label.String = 'Total Intensity of mRNA per Cell';
end

figure;
for jj = 1:size(vorvx,2)
    patch(vorvx{1,jj}(:,1),vorvx{1,jj}(:,2), intens_mRNA(jj) / area_cell(jj));
    set(gca,'YDir','reverse');
    axis equal
    axis(axis_size)
    c = colorbar;
    c.Label.String = 'Total Intensity of mRNA per Unit Area per Cell';
end
