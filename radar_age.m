function [radar] = radar_age(file, cores, Ndraw)

% Conversion to depth
[radar] = radar_depth(file, cores);

% Find the mean response with depth in the radar data attributes across a
% given horizontal resolution (in meters)
horz_res = 25;
[radar] = radar_stack(radar, horz_res);

% Stationarize the radar response by differencing traces with a smoothing 
% spline
s = zeros(size(radar.data_stack));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth(:,i), radar.data_stack(:,i), 0.95, radar.depth(:,i));
end
radar_stat = radar.data_stack - s;

% Remove linear trend in variance (attentuation with depth) and convert to
% standardized values (z-score statistics)
radar_Z = zeros(size(radar_stat));
for i = 1:size(radar_stat, 2)
    data = radar_stat(:,i);
    % Frame length to define local variance
    half_frame = round(0.5*length(data)/5); 
    var0 = movvar(data, 2*half_frame, 'EndPoints', 'discard');
    x = (half_frame:length(data)-half_frame)';
    % Linear trend in variance
    EQ = polyfit(x, var0, 1);
    x_mod = (1:length(data))';
    mod = polyval(EQ, x_mod);
    % Standardize variance-corrected data
    radar_Z(:,i) = data./sqrt(abs(mod));
end

% Define the vertical resolution of the core data
core_res = 0.02;

% Define the cutoff depth for radar traces and find index of crossover
% depth
cutoff = 25;
depth_bott = floor(min([min(radar.depth(end,:)) cutoff]));

% Trim radar traces to cutoff depth and interpolate data to vertical scale
% of the firn cores
radarZ_interp = zeros(depth_bott/core_res+1, size(radar.data_stack, 2));
for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:core_res:radar.depth(end,i));
    radarZ_i = interp1(radar.depth(:,i), radar_Z(:,i), depth_interp, 'pchip');
    radarZ_interp(:,i) = radarZ_i(1:size(radarZ_interp, 1));
end

% If manual layer picks are present, transform them to same depth and
% vertical scale as the interpolated radar data
if isfield(radar, 'man_layers')
    man_interp = zeros(depth_bott/core_res+1, size(radar.data_stack, 2));
    for i = 1:size(radar.data_stack, 2)
        depth_interp = (0:core_res:radar.depth(end,i));
        man_i = interp1(radar.depth(:,i), radar.man_layers(:,i), depth_interp, 'nearest');
        man_interp(:,i) = man_i(1:size(man_interp, 1));
    end
    radar.man_layers = man_interp;
end

% Assign structure output depth to interpolated depths
radar.depth = (0:core_res:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~20 m)
radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);


clearvars -except radar horz_res core_res


% B = ones(5)/5^2;
% data_smooth = conv2(radar.data_smooth, B, 'same');

%% Find depth, width, and prominence of peaks for each radar trace

% Preallocate arrays for various components
peaks_raw = zeros(size(radar.data_smooth));
peak_width = zeros(size(radar.data_smooth));
Proms = cell(1, size(radar.data_smooth, 2));
widths = cell(1,size(radar.data_smooth, 2));
depths = cell(1, size(radar.data_smooth, 2));
depth_idx = cell(1, size(radar.data_smooth, 2));
% troughs = zeros(size(radar.data_smooth));

for i = 1:size(radar.data_smooth, 2)
    data_i = radar.data_smooth(:,i);
%     data_i = data_smooth(:,i);
    
    % Prominence threshold for peaks
    minProm = 0.50;
    % Min distance between peaks (in meters)
    minDist = 0.10;
    
    % Find peak statistics in each trace based on criteria
    [~, peaks_idx_i, widths_i, Prom_i] = findpeaks(data_i, ...
        'MinPeakProminence', minProm, ...
        'MinPeakDistance', minDist/core_res, 'WidthReference', 'halfheight');

    % Add peak prominence and width values to relevent matrices
    peaks_raw(peaks_idx_i,i) = Prom_i;
    peak_width(peaks_idx_i,i) = widths_i;
    
    % Add values to relevent cell arrays
    Proms{i} = Prom_i;
    widths{i} = widths_i;
    depths{i} = radar.depth(peaks_idx_i);
    depth_idx{i} = peaks_idx_i;
    
    
%     [~, trough_idx_i, ~, P_trough_i] = findpeaks(-data_i, ...
%         'MinPeakProminence', minProm, 'MinPeakDistance', minDist/core_res);
%     troughs(trough_idx_i,i) = -P_trough_i;
    
    
end

%%

% Define size of error (in data bins) for radar
% err_bin = round(minDist/core_res);

[~, layers] = find_layers(peaks_raw, peak_width, core_res, horz_res);




%%

% Preallocate arrays for the matrix indices of members of each layer
layers_idx = cell(1,length(layers));
% layer_var = cell(1,length(layers));
peaks = zeros(size(peaks_raw));

% For loop to coerce layers to have one row position for each trace
for i = 1:length(layers_idx)
    
    % Find matrix indices of all members of ith layer
    layer_i = layers{i};
    
    % Find row and col indices of members of ith layer
    [row, col] = ind2sub(size(radar.data_smooth), layer_i);
    mag = peaks_raw(layer_i);
    

    col_interp = min(col):max(col);
    row_interp = round(fnval(csaps(col, row), col_interp));
    row_interp(row_interp < 1) = 1;
    row_interp(row_interp > size(peaks,1)) = size(peaks,1);
    mag_interp = csaps(col, mag, 1/length(col_interp), col_interp);
%     row_var = interp1(col, movvar(row, round(1000/horz_res)), col_interp);
    
    layer_interp = sub2ind(size(peaks), row_interp, col_interp);
    peaks(layer_interp) = mag_interp;
    layers_idx{i} = layer_interp';
%     layer_var = var_interp;
    
%     % If multiple rows exist for the same column, take the
%     % squared prominence-weighted mean of the rows
%     if length(col) > length(unique(col))
%         
%         % Create matrix of peaks with only layer_i members
%         layer_mat = zeros(size(radar.data_smooth));
%         layer_mat(layer_i) = peaks_raw(layer_i);
%         
%         % Find all traces with multiple peaks within layer_i
%         multi_idx = sum(logical(layer_mat))>1;
%         col_nums = 1:size(layer_mat, 2);
%         for k = col_nums(multi_idx)
%             k_idx = find(col==k);
%             k_peaks = layer_mat(row(k_idx),k);
%             
%             % Create squared prominence weighting matrix
%             k_sum = sum(k_peaks.^2);
%             k_row = round(sum((k_peaks.^2/k_sum).*row(k_idx)));
%             k_peak = sum((k_peaks.^2/k_sum).*k_peaks);
% 
%             % Assign squared prominence weighted mean position to layer
%             % trace
%             k_col = zeros(size(layer_mat, 1), 1);
%             k_col(k_row) = k_peak;
%             layer_mat(:,k) = k_col;
%         end
%         peaks_mat = find(layer_mat);
%         peaks_val = layer_mat(peaks_mat);
%         [row, col] = ind2sub(size(radar.data_smooth), peaks_mat);
%         peaks(peaks_mat) = layer_mat(peaks_mat);
%     
%     else
%         % Assign weighted peak position and value to preallocated peaks
%         % matrix
%         peaks(layer_i) = peaks_raw(layer_i);
%     end
    
end

%%

% Calculate global and individual layer reliability for each column in
% radar data
[RMSE_globe, RSE_layer, depth_slope] = REL_score(peaks, layers_idx);


% Calculate continuous layer distances for each layer (accounting for 
% lateral size of stacked radar trace bins)
% layers_dist = cellfun(@(x) numel(x)*horz_res, layers_idx);
layers_dist = cellfun(@(x) numel(x)*horz_res, layers_idx);

% Map layer prominence-distance values to the location within the radar
% matrix of the ith layer
layer_peaks = zeros(size(peaks));
for i = 1:length(layers_idx)
    layer_peaks(layers_idx{i}) = peaks(layers_idx{i}).*layers_dist(i);
end


% Output layer arrays to radar structure
radar.peaks = peaks;
radar.layers = layers_idx;
radar.layer_vals = layer_peaks;

%% Assign layer likelihood scores and estimate age-depth scales

% Define surface age and the year associated with the first pick of the 
% algorithm
age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);

ages = zeros([size(radar.data_smooth) Ndraw]);
radar.likelihood = zeros(size(radar.data_smooth));
for i = 1:size(layer_peaks, 2)
    
%     P_50 = 2*1000;
%     P_50 = 1000*mean(std(radar.data_smooth));
%     P_50 = 1000*(quantile(radar.data_smooth(:,i), 0.95) - ...
%         quantile(radar.data_smooth(:,i), 0.05));
    P_50 = 1*2500*median(Proms{i});
%     P_50 = median(Proms{i})*mean(cellfun(@length, layers_idx));
    
    Po = 0.001;
    K = 1;
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    
    peaks_i = layer_peaks(:,i);
    %     peaks_i = peaks(:,i).*layer_peaks(:,i);

    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar.depth(peaks_idx);
    
    % Probability of peak representing a year based on a logistic function
    % with rate (r) calculated above
    likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i)); % times the prominence at that location?
    radar.likelihood(peaks_idx,i) = likelihood;
    
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end
    
    for j = 1:Ndraw
        depths_j = [0; depths_i(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        ages(:,i,j) = interp1(depths_j, yrs_j, radar.depth, 'linear', 'extrap');
%         try
%             ages(:,i,j) = interp1(depths_j, yrs_j, radar.depth, 'linear', 'extrap');
%         catch ME
%             sprintf('Error age interpolation for trace %u, trial %u. Filling with NAN', i, j);
%             ages(:,i,j) = nan;
%         end
    end
end

radar.ages = ages;

end
