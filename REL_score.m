% Function that looks at the similiarity in relative row position between
% different picked layers (a measure of how effectively/consistently the
% method picked layers in the given data)

function [RMSE, s_matrix] = REL_score(peaks, layers, horz_res)

% Add horizontal surface layer to layer groups
layers{end+1} = sub2ind(size(peaks), ...
    ones(size(peaks,2),1), (1:size(peaks,2))');

% Find the row and column locations of all members of all layers
[~, c_all] = cellfun(@(x) ind2sub(size(peaks), x), ...
    layers, 'UniformOutput', false);

% Define length of data chunk segments over which to estimate layer
% deviation from segment mean as the mean length of picked layers
% seg_length = round(mean(cellfun(@length, layers_idx)));
seg_length = round(500/horz_res);

% Preallocate arrays
depth_slope = zeros(2, size(peaks,2));
% depth_cells = cell(1, size(peaks,2));
RMSE = zeros(1, size(peaks,2));


col_idx = [1 seg_length/2];

% Define column range for data chunk
cols = col_idx(1):col_idx(2);

% Find subset of picked layers which cover the full span of the data
% chunk, and find the row/col subscripts for those subsetted layers
layers_set = layers(cellfun(@(x) min(x) <= col_idx(1) & ...
    max(x) >= col_idx(2), c_all));
[r_set, c_set] = cellfun(@(x) ind2sub(size(peaks), x), ...
    layers_set, 'UniformOutput', false);

% Extract row positions of layer subset that lies within data chunk
% range
rows = cell(1, length(c_set));
for j = 1:length(rows)
    r_idx = c_set{j} >= col_idx(1) & c_set{j} <= col_idx(2);
    rows{j} = r_set{j}(r_idx);
end


p = cellfun(@(y) polyfit(cols, y', 1), rows, 'UniformOutput', false);
slope = cellfun(@(x) x(1), p);
row_depth = cellfun(@mean, rows);
[p, stats] = robustfit(row_depth, slope);
depth_slope(:,cols) = repmat(p, 1, length(cols));
RMSE(cols) = repmat(stats.robust_s, 1, length(cols));





col_i = (seg_length/2)+1:size(peaks,2)-(seg_length/2);

for i = col_i
    cols = i-(seg_length/2):i+(seg_length/2);
    
    % Find subset of picked layers which cover the full span of the data
    % chunk, and find the row/col subscripts for those subsetted layers
    layers_set = layers(cellfun(@(x) min(x) <= cols(1) & ...
        max(x) >= cols(end), c_all));
    [r_set, c_set] = cellfun(@(x) ind2sub(size(peaks), x), ...
        layers_set, 'UniformOutput', false);
    
    % Extract row positions of layer subset that lies within data chunk
    % range
    rows = cell(1, length(c_set));
    for j = 1:length(rows)
        r_idx = c_set{j} >= cols(1) & c_set{j} <= cols(end);
        rows{j} = r_set{j}(r_idx);
    end
    
    p = cellfun(@(y) polyfit(cols, y', 1), rows, 'UniformOutput', false);
    slope = cellfun(@(x) x(1), p);
    row_depth = cellfun(@mean, rows);
    [depth_slope(:,i), stats] = robustfit(row_depth, slope);
    RMSE(i) = stats.robust_s;
    
end


col_idx = [size(peaks,2)-seg_length/2+1 size(peaks,2)];

% Define column range for data chunk
cols = col_idx(1):col_idx(2);

% Find subset of picked layers which cover the full span of the data
% chunk, and find the row/col subscripts for those subsetted layers
layers_set = layers(cellfun(@(x) min(x) <= col_idx(1) & ...
    max(x) >= col_idx(2), c_all));
[r_set, c_set] = cellfun(@(x) ind2sub(size(peaks), x), ...
    layers_set, 'UniformOutput', false);

% Extract row positions of layer subset that lies within data chunk
% range
rows = cell(1, length(c_set));
for j = 1:length(rows)
    r_idx = c_set{j} >= col_idx(1) & c_set{j} <= col_idx(2);
    rows{j} = r_set{j}(r_idx);
end

p = cellfun(@(y) polyfit(cols, y', 1), rows, 'UniformOutput', false);
slope = cellfun(@(x) x(1), p);
row_depth = cellfun(@mean, rows);
[p, stats] = robustfit(row_depth, slope);
depth_slope(:,cols) = repmat(p, 1, length(cols));
RMSE(cols) = repmat(stats.robust_s, 1, length(cols));


depths = (1:size(peaks,1))';
s_matrix = depths*depth_slope(2,:) + depth_slope(1,:);
% s_matrix = depths.*(depth_slope(2,:)+1);

% % Diagnostic plot
% ystart = 10:25:size(peaks,1);
% xstart = ones(1, length(ystart));
% XY = stream2(ones(size(peaks)), s_matrix, xstart, ystart, 1);
% 
% figure
% imagesc(peaks)
% hold on
% hlines = streamline(XY);
% set(hlines, 'LineWidth', 1.5, 'Color', 'r')
% hold off


end

