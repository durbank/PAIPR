function [radar] = radar_age(file, cores, Ndraw)

% Conversion to depth
[radar] = radar_depth(file, cores);

% Find the mean response with depth in the resampled radar data across a
% given lateral distance 'window' (in this case 10 m)
% radar.data_out = movmean(radar.data_out, round(75/mean(diff(radar.dist))), 2);
[radar] = radar_stack(radar, 25);

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
%     radar_Z(:,i) = zscore(data./sqrt(abs(mod)));
    radar_Z(:,i) = data./sqrt(abs(mod));
end

resolution = 0.02;
cutoff = 25;
depth_bott = floor(min([min(radar.depth(end,:)) cutoff]));

radarZ_interp = zeros(depth_bott/resolution+1, size(radar.data_stack, 2));

for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:resolution:radar.depth(end,i));
    radarZ_i = interp1(radar.depth(:,i), radar_Z(:,i), depth_interp, 'pchip');
    radarZ_interp(:,i) = radarZ_i(1:size(radarZ_interp, 1));
end

% Assign output depth to interpolated depths
radar.depth = (0:resolution:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~20 m)
radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);
% radar.data_smooth = sgolayfilt(movmean(radarZ_interp, 5, 2), 3, 9);

% Year associated with the first pick of the algorithm
age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);

%% Find the depth and prominence of peaks within each smoothed radar trace

% Preallocate arrays for various components
peaks_raw = zeros(size(radar.data_smooth));
peak_width = zeros(size(radar.data_smooth));
Proms = cell(1, size(radar.data_smooth, 2));
widths = cell(1,size(radar.data_smooth, 2));
depths = cell(1, size(radar.data_smooth, 2));
depth_idx = cell(1, size(radar.data_smooth, 2));

troughs = zeros(size(radar.data_smooth));

for i = 1:size(radar.data_smooth, 2)
    data_i = radar.data_smooth(:,i);
    minProm = 0.25;                 % Prominence threshold for peaks
    minDist = 0.08;                 % Min distance between peaks (in meters)
    
    % Find peaks in each trace based on requirements
    [~, peaks_idx_i, widths_i, Prom_i] = findpeaks(data_i, ...
        'MinPeakProminence', minProm, ...
        'MinPeakDistance', minDist/resolution, 'WidthReference', 'halfheight');
%     peaks_i = zeros(length(data_i), 1);
%     peaks_i(peaks_idx_i) = Prom_i;
%     peaks(:,i) = peaks_i;

    % Add peak prominence and width values to relevent matrices
    peaks_raw(peaks_idx_i,i) = Prom_i;
    peak_width(peaks_idx_i,i) = widths_i;
    
    % Add values to relevent cell arrays
    Proms{i} = Prom_i;
    widths{i} = widths_i;
    depths{i} = radar.depth(peaks_idx_i);
    depth_idx{i} = peaks_idx_i;
    
    
    [~, trough_idx_i, ~, P_trough_i] = findpeaks(-data_i, ...
        'MinPeakProminence', minProm, 'MinPeakDistance', minDist/resolution);
    troughs(trough_idx_i,i) = -P_trough_i;
    
    
end

%%

% Define size of ~quasi bin confidence interval
err_bin = round(minDist/resolution);

% Preallocate cell array for layer numbers and initialize values by
% assigning unique layer numbers to each peak in the first trace
Groups = Proms;
Groups{1} = uint32((1:length(Groups{1})))';

% Set value for next unique layer number
new_group = Groups{1}(end) + 1;

% Preallocate matrix for layer numbers and add initialized values for the
% first trace (col=1)
peak_group = zeros(size(peaks_raw));
peak_group(depth_idx{1},1) = Groups{1};

for i = 2:size(peaks_raw, 2)
    
    % Assign column bounds for the ith local search window based on 250 m
    % window
    col_idx = [max([i-round(250/mean(diff(radar.dist))) 1]) i-1];
%     col_idx = [max([i-round(0.5*err_bin) 1]) i-1];
    
    for j = 1:length(Proms{i})
        
        % Determine cell array index for jth peak in ith trace
        j_idx = depth_idx{i}(j);
        
        % Determine index values for the row boundaries of the local search
        % window of peak (i,j) based on the bin error window size and the
        % half-width of peak (i,j)
        row_idx = [max([j_idx - 3*round(err_bin+0.5*widths{i}(j)) 1]) ...
            min([j_idx + 3*round(err_bin+0.5*widths{i}(j)) size(peaks_raw, 1)])];
%         row_idx = [max([j_idx-err_bins 1]) min([j_idx+err_bins size(peaks, 1)])];
        
        % Define local window to search for matching layer numbers
        peaks_local = peaks_raw(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        
        % Define the local group matrix
        group_local = peak_group(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        
        % Find the mean row, nearest col, and magnitude of groups within
        % the local window
        group_list = unique(group_local(group_local>0));
        group_row = zeros(length(group_list), 1);
        group_col = zeros(length(group_list), 1);
        group_numel = zeros(length(group_list), 1);
        group_mag = zeros(length(group_list), 1);
        for k = 1:length(group_list)
            [k_rows, k_cols] = find(group_local==group_list(k));
            group_idx = sub2ind(size(group_local), k_rows, k_cols);
            group_numel(k) = length(group_idx);
%             group_val(k) = median(peaks_local(group_idx));
            group_mag(k) = sum(peaks_local(group_idx));
            
%             weights = (peaks_local(group_idx) + k_cols)/...
%                 sum(peaks_local(group_idx) + k_cols);
%             group_row(k) = sum(weights.*k_rows);
            if group_numel(k) > 3
                p = polyfit(k_cols, k_rows, 1);
                group_row(k) = polyval(p, size(group_local, 2)+1);
            else
                group_row(k) = mean(k_rows);
            end
            group_col(k) = max(k_cols) + 1;
        end
        
        
        
        
        
        
%         % Calculate distances between peak (i,j) and mean group values
%         % within local window based on differences in depth, lateral
%         % distance, and peak prominence
%         w_dist = sqrt((j_idx - (row_idx(1)+group_row-1)).^2 + ...
%             (i - (col_idx(1)+group_col-1)).^2 + ...
%             (Proms{i}(j)-group_mag./group_numel).^2);
%         
%         % Assign distance threshold based on the (i,j) peak half-width and
%         % the error bin size
%         threshold = (0.5*widths{i}(j) + err_bin);
%         
%         % Create logical index of local peak groups within threshold
%         % tolerance, and apply to the group_mag array
%         tol_idx = w_dist <= threshold;
%         group_mag = group_mag(tol_idx);
%         
%         % Scale local group distances within tolerance using the inverse
%         % group magnitude value (stronger peaks have smaller distances than
%         % equivalent weaker peaks) and find peak group with the shortest
%         % distance
%         [~, dist_idx] = min((1./group_mag).*w_dist(tol_idx));
%         
%         if isempty(dist_idx)
%             % If peak (i,j) nearest neighbor has a distance greater than
%             % the threshold, assign a new unqiue layer number to peak (i,j)
%             Groups{i}(j) = new_group;
%             peak_group(j_idx,i) = new_group;
%             new_group = new_group + 1;
%             
% 
%         else
%             % Assign peak (i,j) to the nearest neighbor group
%             group_j = uint32(group_list(dist_idx));
%             Groups{i}(j) = group_j;
%             peak_group(j_idx,i) = group_j;



        
        % Calculate distances between peak (i,j) and mean group values
        % within local window based on differences in depth, lateral
        % distance, and peak prominence
%         w_dist = sqrt((j_idx - (row_idx(1)+group_row-1)).^2 + ...
%             (i - (col_idx(1)+group_col-1)).^2 + ...
%             (Proms{i}(j)-group_val./group_numel).^2);
        w_dist = (1./group_mag).*sqrt((j_idx - (row_idx(1)+group_row-1)).^2 + ...
            (i - (col_idx(1)+group_col-1)).^2 + ...
            (Proms{i}(j)-group_mag./group_numel).^2);
        
        % Select the nearest group neighbor to peak (i,j)
        [min_dist, dist_idx] = min(w_dist);
        
        % Set distance threshold based on peaks 50 m laterally apart and
        % the error bin size, scaled by the (i,j) peak prominence
        threshold = (0.5*widths{i}(j) + err_bin);
            
        if min_dist.*group_mag(dist_idx) <= threshold
            % Assign peak (i,j) to the nearest neighbor group
            group_j = uint32(group_list(dist_idx));
            Groups{i}(j) = group_j;
            peak_group(j_idx,i) = group_j;

        else
            % If peak (i,j) nearest neighbor has a distance greater than
            % the threshold, assign a new unqiue layer number to peak (i,j)
            Groups{i}(j) = new_group;
            peak_group(j_idx,i) = new_group;
            new_group = new_group + 1;
            
            
            
            
        end
    end
end





% % Preallocate arrays for the matrix indices of members of each layer
% layers_idx = cell(1,new_group-1);
% 
% for i = 1:length(layers_idx)
%     
%     % Find matrix indices of all members of ith layer, and assign to
%     % preallocated cell array
%     layers_idx{i} = find(peak_group==i);
% end
% 
% layers_idx = layers_idx(cellfun(@(x) ...
%     length(x) > round(100/mean(diff(radar.dist))), layers_idx));







% Preallocate arrays for the matrix indices of members of each layer
layers_idx = cell(1,new_group-1);
peaks = zeros(size(peaks_raw));
for i = 1:length(layers_idx)
    
    % Find matrix indices of all members of ith layer, and assign to
    % preallocated cell array
    layer_i = find(peak_group==i);
    
    % Find row and col indices of members of ith layer
    [row, col] = ind2sub(size(radar.data_smooth), layer_i);
    
    % If multiple rows exist for the same column, take the
    % squared prominence-weighted mean of the rows
    if length(col) > length(unique(col))
        layer_mat = zeros(size(radar.data_smooth));
        layer_mat(layer_i) = peaks_raw(layer_i);
        multi_idx = sum(logical(layer_mat))>1;
        col_nums = 1:size(layer_mat, 2);
        for k = col_nums(multi_idx)
            k_idx = find(col==k);
            k_peaks = layer_mat(row(k_idx),k);
            k_sum = sum(k_peaks.^2);
            k_row = round(sum((k_peaks.^2/k_sum).*row(k_idx)));
            k_peak = sum((k_peaks.^2/k_sum).*k_peaks);
%             k_sum = sum(k_peaks);
%             k_row = round(sum((k_peaks/k_sum).*row(k_idx)));
%             k_peak = sum((k_peaks/k_sum).*k_peaks);
            k_col = zeros(size(layer_mat, 1), 1);
            k_col(k_row) = k_peak;
            layer_mat(:,k) = k_col;
        end
        peaks_mat = find(layer_mat);
        [row, col] = ind2sub(size(radar.data_smooth), peaks_mat);
        peaks(peaks_mat) = layer_mat(peaks_mat);
    
    else
        peaks(layer_i) = peaks_raw(layer_i);
    end

%     % Smooth layer i using a moving average of row indices
%     row_mean = round(movmean(row, 10));
%     layers_test{i} = sub2ind(size(radar.data_smooth), row_mean, col);
    layers_idx{i} = sub2ind(size(radar.data_smooth), row, col);
    
%     row_mean = round(movmean(row, 10));
%     
%     % If multiple rows exist for the same column, select the nearest
%     % neighbor compared to the rest of the layer, and remove others
%     if length(col) > length(unique(col))
%         layer_mat = zeros(size(radar.data_smooth));
%         layer_mat(layers_idx{i}) = peaks(layers_idx{i});
%         multi_idx = sum(logical(layer_mat))>1;
%         col_nums = 1:size(layer_mat, 2);
%         for k = col_nums(multi_idx)
%             k_col = layer_mat(:,k);
%             k_idx = find(k_col);
%             P_bounds = [max([k-10 1]) min([k+10 size(layer_mat, 2)])];
%             P_k = layer_mat(:,P_bounds(1):P_bounds(2));
%             P_k = median(P_k(P_k>0));
%             k_dist = sqrt((round(mean(row_mean(col==k))) - k_idx).^2 + ...
%                 (P_k - k_col(k_idx)).^2);
%             [~,min_idx] = min(k_dist);
%             k_col = zeros(size(layer_mat, 1), 1);
%             k_col(k_idx(min_idx)) = layer_mat(k_idx(min_idx),k);
%             layer_mat(:,k) = k_col;
%         end
%         layers_idx{i} = find(layer_mat);
%     end
    
end


% Integrate peak magnitudes across ith layer to obtain layer
% prominence-distance value (accounting for lateral size of stacked
% radar trace bins)
% layers_val = cellfun(@(x) sum(peaks(x))*mean(diff(radar.dist)), layers);
layers_dist = cellfun(@(x) numel(x)*median(diff(radar.dist)), layers_idx);

% Map layer prominence-distance values to the location within the radar
% matrix of the ith layer
layer_peaks = zeros(size(peaks));
for i = 1:length(layers_idx)
    layer_peaks(layers_idx{i}) = peaks(layers_idx{i}).*layers_dist(i);
end

layers_idx = layers_idx(cellfun(@(x) ...
    length(x) > round(100/mean(diff(radar.dist))), layers_idx));

%%

% Output layer arrays to radar structure
radar.peaks = peaks;
radar.layers = layers_idx;
radar.layer_vals = layer_peaks;



ages = zeros([size(radar.data_smooth) Ndraw]);
radar.likelihood = zeros(size(radar.data_smooth));
for i = 1:size(layer_peaks, 2)
    
%     P_50 = 2*1000;
%     P_50 = 1000*mean(std(radar.data_smooth));
%     P_50 = 1000*(quantile(radar.data_smooth(:,i), 0.95) - ...
%         quantile(radar.data_smooth(:,i), 0.05));
    P_50 = 1*500*median(Proms{i});
    
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
    end
end

radar.ages = ages;

end
