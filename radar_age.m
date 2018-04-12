function [radar] = radar_age(file, cores, Ndraw)

% Conversion to depth
[radar] = radar_depth(file, cores);

% Find the mean response with depth in the resampled radar data across a
% given lateral distance 'window' (in this case ~100 m)
[radar] = radar_stack(radar, 30);

% Stationarize the radar response using a smoothing spline
s = zeros(size(radar.data_stack));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth(:,i), radar.data_stack(:,i), 0.95, radar.depth(:,i));
end
radar_stat = radar.data_stack - s;

% Remove linear trend in variance (attentuation with depth) and convert to
% z-score statistics
radar_Z = zeros(size(radar_stat));
for i = 1:size(radar_stat, 2)
    data = radar_stat(:,i);
    half_frame = round(0.5*length(data)/5);
    var0 = movvar(data, 2*half_frame, 'EndPoints', 'discard');
    x = (half_frame:length(data)-half_frame)';
    EQ = polyfit(x, var0, 1);
    x_mod = (1:length(data))';
    mod = EQ(1)*x_mod+EQ(2);
    radar_Z(:,i) = zscore(data./sqrt(abs(mod)));
end

resolution = 0.02;
depth_bott = floor(min([min(radar.depth(end,:)) 30]));

radarZ_interp = zeros(depth_bott/resolution+1, size(radar.data_stack, 2));

for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:resolution:radar.depth(end,i));
    radarZ_i = interp1(radar.depth(:,i), radar_Z(:,i), depth_interp, 'pchip');
    radarZ_interp(:,i) = radarZ_i(1:size(radarZ_interp, 1));
end

radar.depth = (0:resolution:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~20 m)
radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);

% Year associated with the first pick of the algorithm
% age_top = round(radar.collect_date);
age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);

%% Find the depth and prominence of peaks within each smoothed radar trace

% Preallocate arrays for various components
peaks = zeros(size(radar.data_smooth));
peak_width = zeros(size(radar.data_smooth));
Proms = cell(1, size(radar.data_smooth, 2));
widths = cell(1,size(radar.data_smooth, 2));
depths = cell(1, size(radar.data_smooth, 2));
depth_idx = cell(1, size(radar.data_smooth, 2));

for i = 1:size(radar.data_smooth, 2)
    data_i = radar.data_smooth(:,i);
    minProm = 0.05;                 % Prominence threshold for peaks
%     minProm = iqr(data_i);
%     minDist = max([0 min(dist_Ex)-3*dist_STD]);
    minDist = 0.12;                 % Min distance between peaks (in meters)
    
    % Find peaks in each trace based on requirements
    [~, peaks_idx_i, widths_i, Prom_i] = findpeaks(data_i, ...
        'MinPeakProminence', minProm, ...
        'MinPeakDistance', minDist/resolution);
%     peaks_i = zeros(length(data_i), 1);
%     peaks_i(peaks_idx_i) = Prom_i;
%     peaks(:,i) = peaks_i;

    % Add peak prominence and width values to relevent matrices
    peaks(peaks_idx_i,i) = Prom_i;
    peak_width(peaks_idx_i,i) = widths_i;
    
    % Add values to relevent cell arrays
    Proms{i} = Prom_i;
    widths{i} = widths_i;
    depths{i} = radar.depth(peaks_idx_i);
    depth_idx{i} = peaks_idx_i;
end

%%

% Define size of ~quasi bin confidence interval
err_bin = minDist/resolution;
err_bin = 4;

% Preallocate cell array for layer numbers and initialize values by
% assigning unique layer numbers to each peak in the first trace
Groups = Proms;
Groups{1} = uint32((1:length(Groups{1})))';

% Set value for next unique layer number
new_group = Groups{1}(end) + 1;

% Preallocate matrix for layer numbers and add initialized values for the
% first trace (col=1)
peak_group = zeros(size(peaks));
peak_group(depth_idx{1},1) = Groups{1};

for i = 2:size(peaks, 2)
    
    % Assign column bounds for the ith local search window based on 100 m
    % window
    col_idx = [max([i-round(100/mean(diff(radar.dist))) 1]) i-1];
%     col_idx = [max([i-round(0.5*err_bin) 1]) i-1];
    
    for j = 1:length(Proms{i})
        
        % Determine cell array index for jth peak in ith trace
        j_idx = depth_idx{i}(j);
        
        % Determine index values for the row boundaries of the local search
        % window of peak (i,j) based on the bin error window size and the
        % half-width of peak (i,j)
        row_idx = [max([j_idx-round(0.5*(err_bin+widths{i}(j))) 1]) ...
            min([j_idx+round(0.5*(err_bin+widths{i}(j))) size(peaks, 1)])];
        
        % Define local window to search for matching layer numbers
        peaks_local = peaks(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        
        % Find the row, col, and index values for peaks within the local
        % search window
        [local_row, local_col] = find(peaks_local);
        local_idx = sub2ind(size(peaks_local), local_row, local_col);
        
        % Calculate distances between peak (i,j) and peaks within the local
        % search window based on differences in depth, lateral distance,
        % and peak prominence magnitudes
        w_dist = sqrt((j_idx - (row_idx(1)+local_row-1)).^2 + ...
            (i - (col_idx(1)+local_col-1)).^2 + ...
            (Proms{i}(j) - peaks_local(local_idx)).^2);
        
        % Select the nearest neighbor to peak (i,j)
        [~, dist_idx] = min(w_dist);
        
        % (I may add a tolerance in the future)
        if ~isempty(dist_idx)     % if j_min <= bin_res/2
            
            % If present, assign the neighest neighbor layer number to the 
            % peak (i,j) layer number in both cell array and matrix
            group_j = uint32(peak_group(row_idx(1)+local_row(dist_idx)-1,...
                col_idx(1)+local_col(dist_idx)-1));
            Groups{i}(j) = group_j;
            peak_group(j_idx,i) = group_j;
        else
            
            % If peak (i,j) does not have a nearest neighbor, assign a new
            % unqiue layer number to peak (i,j)
            Groups{i}(j) = new_group;
            peak_group(j_idx,i) = new_group;
            new_group = new_group + 1;
        end
    end
end

% Preallocate arrays for the matrix indices of members of each layer and
% integrated prominence-distance mapped to the radar matrix
layers_idx = cell(1,new_group-1);
layer_peaks = zeros(size(peaks));
for i = 1:length(layers_idx)
    
    % Find matrix indices of all members of ith layer, and assign to
    % preallocated cell array
    layers_idx{i} = find(peak_group==i);
    
    % Find row and col indices of members of ith layer
    [~, col] = ind2sub(size(radar.data_smooth), layers_idx{i});
    
    % If multiple rows exist for the same column, select the strongest
    % peak among them and remove others
    if length(col) > length(unique(col))
        layer_mat = zeros(size(radar.data_smooth));
        layer_mat(layers_idx{i}) = peaks(layers_idx{i});
        multi_idx = sum(logical(layer_mat))>1;
        col_nums = 1:size(layer_mat, 2);
        for k = col_nums(multi_idx)
            k_col = zeros(size(layer_mat, 1), 1);
            [~,k_max] = max(layer_mat(:,k));
            k_col(k_max) = layer_mat(k_max,k);
            layer_mat(:,k) = k_col;
        end
        layers_idx{i} = find(layer_mat);
    end

    % Integrate peak magnitudes across ith layer to obtain layer
    % prominence-distance value (accounting for lateral size of stacked
    % radar trace bins)
%     layers_val(i) = sum(peaks(layers_idx{i}));
    layers_val = sum(peaks(layers_idx{i}))*mean(diff(radar.dist));
    
    % Map layer prominence-distance values to the location within the radar
    % matrix of the ith layer
    layer_peaks(layers_idx{i}) = layers_val;
end

% Output layer arrays to radar structure
radar.layers = layers_idx;
radar.layer_vals = layer_peaks;

% P_50 = 2*1000;
% P_50 = 1000*2*mean(iqr(radar.data_smooth));
P_50 = 1000*mean(std(radar.data_smooth));

Po = 0.001;
K = 1;
r = log((K*Po/0.50-Po)/(K-Po))/-P_50;

ages = zeros([size(radar.data_smooth) Ndraw]);
radar.likelihood = zeros(size(radar.data_smooth));
for i = 1:size(layer_peaks, 2)
    
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
    end
end

radar.age = ages;

end
