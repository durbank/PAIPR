% Function that looks at the similiarity in relative row position between
% different picked layers (a measure of how effectively/consistently the
% method picked layers in the given data)

function [RMSE_globe, depth_slope] = REL_score2(peaks, layers_idx)

% Add horizontal surface layer to layer groups
layers_idx{end+1} = sub2ind(size(peaks), ...
    ones(size(peaks,2),1), (1:size(peaks,2))');

% Find the row and column locations of all members of all layers
[~, c_all] = cellfun(@(x) ind2sub(size(peaks), x), ...
    layers_idx, 'UniformOutput', false);

% Define length of data chunk segments over which to estimate layer
% deviation from segment mean as the mean length of picked layers
seg_length = round(mean(cellfun(@length, layers_idx)));

% Find the column indices of segment endpoints
col_idx = [1:seg_length:size(peaks,2) size(peaks,2)];

% Preallocate arrays
depth_slope = zeros(1, size(peaks,2));
% depth_cells = cell(1, size(peaks,2));
RMSE_globe = zeros(1, size(peaks,2));

% Initialize diagnostic figure
figure
hold on
set(gca, 'Ydir', 'reverse')

% Calculate data structure within each vertical data chunk
for i = 1:length(col_idx)-1
    
    % Define column range for data chunk
    cols = col_idx(i):col_idx(i+1);
    
    % Find subset of picked layers which cover the full span of the data
    % chunk, and find the row/col subscripts for those subsetted layers
    layers_set = layers_idx(cellfun(@(x) min(x) <= col_idx(i) & ...
        max(x) >= col_idx(i+1), c_all));
    [r_set, c_set] = cellfun(@(x) ind2sub(size(peaks), x), ...
        layers_set, 'UniformOutput', false);
    
    % Extract row positions of layer subset that lies within data chunk
    % range
    rows = cell(1, length(c_set));
    for j = 1:length(rows)
        r_idx = c_set{j} >= col_idx(i) & c_set{j} <= col_idx(i+1);
        rows{j} = r_set{j}(r_idx);
    end
    
    % Find the mean depth of clipped layer subsets, and sort from surface
    % to the deepest layer
    [row_depth, sort_idx] = sort(cellfun(@mean, rows'));
    
    % Find the deviation of each row position from the layer mean within
    % each clipped layer subset and order by mean depth
    row_mat = (cell2mat(cellfun(@(x) x - mean(x), rows, 'UniformOutput', false)))';
    row_dev = row_mat(sort_idx,:);
    
    % Calculate the best-fit deviation of row position (relative to 
    % respective layer mean row position) with depth for each column within
    % data chunk, and estimate the global RMSE of our layer model for each
    % column in data chunk
    for j = 1:length(cols)
        [p, stats] = robustfit(row_depth, row_dev(:,j));
        depth_slope(col_idx(i)+j-1) = p(2);
        RMSE_globe(col_idx(i)+j-1) = stats.robust_s;
    end
    
%     depth_cells(col_idx(i):col_idx(i+1)) = {row_depth};
    
    % Diagnostic plot
    mdl = (depth_slope(col_idx(i):col_idx(i+1))+1).*row_depth;
    for j = 1:length(rows)
        plot(cols, rows{j}, 'LineWidth', 2)
        plot(cols, mdl(j,:), 'k')
    end
    
end

hold off

%% Individual reliability scores

% [r_all, c_all] = cellfun(@(x) ind2sub(size(peaks), x), layers_idx, 'UniformOutput', false);
% row_depth = cellfun(@mean, r_all');
% 
% RSE_layer = cell(1, length(layers_idx));
% for i = 1:length(RSE_layer)
%     mdl = (depth_slope(c_all{i}(1):c_all{i}(end))+1).*row_depth(i);
%     RSE_layer{i} = sqrt((r_all{i} - mdl').^2);
% end
    




end