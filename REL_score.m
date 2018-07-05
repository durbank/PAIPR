% Function that looks at the similiarity in relative row position between
% different picked layers (a measure of how effectively/consistently the
% method picked layers in the given data)

function [reliability, RMSE, s_matrix] = REL_score(peaks, layers, horz_res)

% Add horizontal surface layer to layer groups
layers{end+1} = sub2ind(size(peaks), ...
    ones(size(peaks,2),1), (1:size(peaks,2))');

% Find the row and column subscripts of all members for all layers
[~, c_all] = cellfun(@(x) ind2sub(size(peaks), x), ...
    layers, 'UniformOutput', false);

% Define length of data chunk segments over which to estimate layer
% layer slopes
seg_length = round(500/horz_res);
% seg_length = round(mean(cellfun(@length, layers_idx)));

% Preallocate arrays
depth_slope = zeros(2, size(peaks,2));
RMSE = zeros(1, size(peaks,2));
R = zeros(1, size(peaks,2));

%% Calculations for first chunk of data (non-iterative)

% Define data chunk endpoints
col_idx = [1 seg_length/2];

% Define column range for data chunk
cols = col_idx(1):col_idx(2);

% Find subset of picked layers which cover the full span of the data
% chunk, and find the row/col subscripts for that layer subset
layers_set = layers(cellfun(@(x) min(x) <= col_idx(1) & ...
    max(x) >= col_idx(2), c_all));
[r_set, c_set] = cellfun(@(x) ind2sub(size(peaks), x), ...
    layers_set, 'UniformOutput', false);

% Extract row indices of layer subset that lies within data chunk range
rows = cell(1, length(c_set));
for j = 1:length(rows)
    r_idx = c_set{j} >= col_idx(1) & c_set{j} <= col_idx(2);
    rows{j} = r_set{j}(r_idx);
end

% Calculate the slope of each layer within the subset for the current data
% chunk using 1st order polynomial fitting
p = cellfun(@(y) polyfit(cols, y', 1), rows, 'UniformOutput', false);
slope = cellfun(@(x) x(1), p);

% Calculate mean depth for each layer subset
row_depth = cellfun(@mean, rows);

% Calculate rate of change in layer slope with depth using linear
% regression robust against outliers
[p, stats] = robustfit(row_depth, slope);
depth_slope(:,cols) = repmat(p, 1, length(cols));

% Output the root-mean-squared error for the rate of change in layer slope
% with depth
RMSE(cols) = repmat(stats.s, 1, length(cols));
R(cols) = repmat(stats.coeffcorr(2,1), 1, length(cols));

%% Iterative calculations for each trace of radargram data

% Define center points at which to calculate slopes
col_i = (seg_length/2)+1:size(peaks,2)-(seg_length/2);

% Iterate for each centerpoint of data
for i = col_i
    
    % Define window over which to calculate layer slopes
    cols = i-(seg_length/2):i+(seg_length/2);
    
    % Find subset of picked layers which cover the full span of the data
    % chunk, and find the row/col subscripts for that layer subset
    layers_set = layers(cellfun(@(x) min(x) <= cols(1) & ...
        max(x) >= cols(end), c_all));
    [r_set, c_set] = cellfun(@(x) ind2sub(size(peaks), x), ...
        layers_set, 'UniformOutput', false);
    
    % Extract row indices of layer subset that lies within data chunk range
    rows = cell(1, length(c_set));
    for j = 1:length(rows)
        r_idx = c_set{j} >= cols(1) & c_set{j} <= cols(end);
        rows{j} = r_set{j}(r_idx);
    end
    
    % Calculate the slope of each layer within the subset for the current 
    % data chunk using 1st order polynomial fitting
    p = cellfun(@(y) polyfit(cols, y', 1), rows, 'UniformOutput', false);
    slope = cellfun(@(x) x(1), p);
    
    % Calculate mean depth for each layer subset
    row_depth = cellfun(@mean, rows);
    
    % Calculate rate of change in layer slope with depth using linear
    % regression robust against outliers
    [depth_slope(:,i), stats] = robustfit(row_depth, slope);
    
    % Output the root-mean-squared error for the rate of change in layer 
    % slope with depth
    RMSE(i) = stats.s;
    R(i) = stats.coeffcorr(2,1);
end

%% Calculations for final chunk of data (non-iterative)

% Define data chunk endpoints
col_idx = [size(peaks,2)-seg_length/2+1 size(peaks,2)];

% Define column range for data chunk
cols = col_idx(1):col_idx(2);

% Find subset of picked layers which cover the full span of the data
% chunk, and find the row/col subscripts for that layer subset
layers_set = layers(cellfun(@(x) min(x) <= col_idx(1) & ...
    max(x) >= col_idx(2), c_all));
[r_set, c_set] = cellfun(@(x) ind2sub(size(peaks), x), ...
    layers_set, 'UniformOutput', false);

% Extract row indices of layer subset that lies within data chunk range
rows = cell(1, length(c_set));
for j = 1:length(rows)
    r_idx = c_set{j} >= col_idx(1) & c_set{j} <= col_idx(2);
    rows{j} = r_set{j}(r_idx);
end

% Calculate the slope of each layer within the subset for the current data
% chunk using 1st order polynomial fitting
p = cellfun(@(y) polyfit(cols, y', 1), rows, 'UniformOutput', false);
slope = cellfun(@(x) x(1), p);

% Calculate mean depth for each layer subset
row_depth = cellfun(@mean, rows);

% Calculate rate of change in layer slope with depth using linear
% regression robust against outliers
[p, stats] = robustfit(row_depth, slope);
depth_slope(:,cols) = repmat(p, 1, length(cols));

% Output the root-mean-squared error for the rate of change in layer slope
% with depth
RMSE(cols) = repmat(stats.s, 1, length(cols));
R(cols) = repmat(stats.coeffcorr(2,1), 1, length(cols));

% Generate matrix of estimated/interpolated layer slopes for each data
% point in radargram
depths = (1:size(peaks,1))';
% s_matrix = depths*depth_slope(2,:) + depth_slope(1,:);
s_matrix = depths*sgolayfilt(depth_slope(2,:), 2, 2*round((length(R)/10)/2)-1);
RMSE = sgolayfilt(RMSE, 2, 2*round((length(R)/10)/2)-1);
reliability = abs(sgolayfilt(R, 2, 2*round((length(R)/10)/2)-1));
% s_matrix = depths.*(depth_slope(2,:)+1);


end

