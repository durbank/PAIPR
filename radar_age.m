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

% Define surface age and the year associated with the first pick of the 
% algorithm
age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);

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
err_bin = round(minDist/core_res);

% Preallocate layer group matrix and cell array
peak_group = zeros(size(peaks_raw));
layers = cell(1, numel(peaks_raw(peaks_raw>0)));

% Matrix of index values (to be used in logical indexing to replace slower
% 'find' function
matrix_idx = reshape(1:numel(peaks_raw), size(peaks_raw));

% Initialize while loop values
peak_pool = peaks_raw;
group_num = 1;
i = 1;
search_new = true;





% While loop that iterates on each accumulation layer
while search_new == true
    
    % Find the maximum peak remaining in pool and assign to group
    [~, peak_max] = max(peak_pool(:));
    peak_group(peak_max) = group_num;
    
    % Initialize values for group search within layer to right of max peak
    peak_n = peak_max;
    search_R = true;
    
    
    while search_R == true
        
        % Find all current layer group members
        layer_i = matrix_idx(peak_group==group_num);
        
        % Get nearest 5 group members
        group_idx = layer_i(end-min([5 length(layer_i)])+1:end);
        
        % Find subscripts of farthest right layer group member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        
        % Define estimated row position for projected nearest neighbor as 
        % the weighted average of nearest 5 layer members and the farthest
        % right layer member (50/50 weighting)
        [row_temp, ~] = ind2sub(size(peaks_raw), group_idx);
        row_n = round((1/2)*(row_n + mean(row_temp)));
        
        % Define estimated peak magnitude and peak width for the projected  
        % nearest neighbor as the mean magnitude and width of the nearest 
        % 5 group members
        mag_n = mean(peaks_raw(group_idx));
        width_n = mean(peak_width(group_idx));
        
        % Define local  search window as 250 m to right of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_n-round(0.50/core_res) 1]) ...
            min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [min([col_n+1 size(peaks_raw, 2)]) ...
            min([col_n+round(250/horz_res) size(peaks_raw, 2)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        % Calculate distance between projected layer position and peaks
        % within the local search window
        dist_n = sqrt((row_n - (row_idx(1)+row_local-1)).^2 + ...
            2*(col_n - (col_idx(1)+col_local-1)).^2 + (mag_n - mag_local).^2);
        
        % Select the nearest neighbor to the estimated layer position
        [min_dist, dist_idx] = min(dist_n);
        
        % Set distance threshold based on peak width and error bin size
        threshold = 1*width_n + err_bin;
        
        if min_dist <= threshold
            % Find index of nearest neighbor (when within tolerance)
            peak_near = sub2ind(size(peaks_raw), ...
                row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
            
            % Assign nearest neighbor to current layer group and remove
            % from the pool of available peaks to search
            peak_group(peak_near) = group_num;
            peak_pool(peak_near) = 0;
            
            % Set newest layer member as the next peak to use in search
            peak_n = peak_near;
            
            % If layer has reached the right edge of radargram, end search
            % to right of layer
            if (col_idx(1)+col_local(min_idx)-1) >= size(peaks_raw, 2)
                search_R = false;
            end
            
        else
            % If no neighbors are within tolerance, end search to right of
            % layer
            search_R = false;
        end
    end
    
    % Search to right of most recent group member for new group members
    while search_R == true % && j >= 4
        
        % Find all current layer group members
        layer_i = matrix_idx(peak_group==group_num);
        
        [row_i, col_i] = ind2sub(size(peaks_raw), layer_i);
        mag_i = peaks_raw(layer_i);
        
        % Get nearest 5 group members
        group_idx = layer_i(end-min([5 length(layer_i)])+1:end);
        
        % Find subscripts of farthest right layer group member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        
        % Define estimated row position for projected nearest neighbor as 
        % the weighted average of nearest 5 layer members and the farthest
        % right layer member (50/50 weighting)
        [row_temp, ~] = ind2sub(size(peaks_raw), group_idx);
        row_n = round((1/2)*(row_n + mean(row_temp)));
        
        % Define estimated peak magnitude and peak width for the projected  
        % nearest neighbor as the mean magnitude and width of the nearest 
        % 5 group members
        mag_n = mean(peaks_raw(group_idx));
        width_n = mean(peak_width(group_idx));
        
        % Define local  search window as 250 m to right of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_n-round(0.50/core_res) 1]) ...
            min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [min([col_n+1 size(peaks_raw, 2)]) ...
            min([col_n+round(250/horz_res) size(peaks_raw, 2)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        data_i = [row_i col_i mag_i];
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1) mag_local];
        
%         dist_n = diag((data_local - repmat(mean(data_i),length(local_idx), 1))...
%             *inv(cov(data_i))*(data_local - repmat(mean(data_i),length(local_idx),1))');
        
        
        % Calculate distance between projected layer position and peaks
        % within the local search window
        dist_n = sqrt((row_n - (row_idx(1)+row_local-1)).^2 + ...
            2*(col_n - (col_idx(1)+col_local-1)).^2 + (mag_n - mag_local).^2);
        
        % Select the nearest neighbor to the estimated layer position
        [min_dist, dist_idx] = min(dist_n);
        
        % Set distance threshold based on peak width and error bin size
        threshold = 1*width_n + err_bin;
        
        % Determine if nearest neighbor is within distance tolerance
        if min_dist <= threshold
            
            % Find index of nearest neighbor (when within tolerance)
            peak_near = sub2ind(size(peaks_raw), ...
                row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
            
            % Assign nearest neighbor to current layer group and remove
            % from the pool of available peaks to search
            peak_group(peak_near) = group_num;
            peak_pool(peak_near) = 0;
            
            % Set newest layer member as the next peak to use in search
            peak_n = peak_near;
            
            % If layer has reached the right edge of radargram, end search
            % to right of layer
            if (col_idx(1)+col_local(min_idx)-1) >= size(peaks_raw, 2)
                search_R = false;
            end
            
        else
            % If no neighbors are within tolerance, end search to right of
            % layer
            search_R = false;
        end
    end
    
    % Initialize values for group search within layer to left of max peak
    peak_n = peak_max;
    search_L = true;
    
    for j = 1:4
        % Find all current layer group members
        layer_i = matrix_idx(peak_group==group_num);
        
        % Get nearest 5 group members
        group_idx = layer_i(1:min([5 length(layer_i)]));
        
        % Find subscripts of farthest left layer member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        
        % Define estimated row position for projected nearest neighbor as
        % the weighted average of nearest 5 layer members and the farthest
        % left layer member (50/50 weighting)
        [row_i, ~] = ind2sub(size(peaks_raw), group_idx);
        row_n = round((1/2)*(row_n + mean(row_i)));
        
        % Define estimated peak magnitude and peak width for the projected
        % nearest neighbor as the mean magnitude and width of the nearest
        % 5 group members
        mag_n = mean(peaks_raw(group_idx));
        width_n = mean(peak_width(group_idx));
        
        % Define local  search window as 250 m to left of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_n-round(0.50/core_res) 1]) ...
            min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col_n-round(250/horz_res)-1]) col_n-1];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        % Calculate distance between projected layer position and peaks
        % within the local search window
        dist_n = sqrt((row_n - (row_idx(1)+row_local-1)).^2 + ...
            2*(col_n - (col_idx(1)+col_local-1)).^2 + (mag_n - mag_local).^2);
        
        % Select the nearest neighbor to the estimated layer position
        [min_dist, dist_idx] = min(dist_n);
        
        % Set distance threshold based on peak width and error bin size
        threshold = 1*width_n + err_bin;
        
        % Determine if nearest neighbor is within distance tolerance
        if min_dist <= threshold
            
            % Find index of nearest neighbor (when within tolerance)
            peak_near = sub2ind(size(peaks_raw), ...
                row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
            
            % Assign nearest neighbor to current layer group and remove
            % from the pool of available peaks to search
            peak_group(peak_near) = group_num;
            peak_pool(peak_near) = 0;
            
            % Set newest layer member as the next peak to use in search
            peak_n = peak_near;
            
        else
            % If no neighbors are within tolerance, end search to left of
            % layer
            search_L = false;
        end
    end
    
    
    
    
    
    % Search to left of most recent group member for new group members
    while search_L == true % && j >= 4
        
        % Find all current layer group members
        layer_i = matrix_idx(peak_group==group_num);
        
        % Get nearest 5 group members
        group_idx = layer_i(1:min([5 length(layer_i)]));
        
        % Find subscripts of farthest left layer member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        
        % Define estimated row position for projected nearest neighbor as 
        % the weighted average of nearest 5 layer members and the farthest
        % left layer member (50/50 weighting)
        [row_i, ~] = ind2sub(size(peaks_raw), group_idx);
        row_n = round((1/2)*(row_n + mean(row_i)));
        
        % Define estimated peak magnitude and peak width for the projected  
        % nearest neighbor as the mean magnitude and width of the nearest 
        % 5 group members
        mag_n = mean(peaks_raw(group_idx));
        width_n = mean(peak_width(group_idx));
        
        % Define local  search window as 250 m to left of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_n-round(0.50/core_res) 1]) ...
            min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col_n-round(250/horz_res)-1]) col_n-1];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        % Calculate distance between projected layer position and peaks
        % within the local search window
        dist_n = sqrt((row_n - (row_idx(1)+row_local-1)).^2 + ...
            2*(col_n - (col_idx(1)+col_local-1)).^2 + (mag_n - mag_local).^2);
        
        % Select the nearest neighbor to the estimated layer position
        [min_dist, dist_idx] = min(dist_n);
        
        % Set distance threshold based on peak width and error bin size
        threshold = 1*width_n + err_bin;
        
        % Determine if nearest neighbor is within distance tolerance
        if min_dist <= threshold
            
            % Find index of nearest neighbor (when within tolerance)
            peak_near = sub2ind(size(peaks_raw), ...
                row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
            
            % Assign nearest neighbor to current layer group and remove
            % from the pool of available peaks to search
            peak_group(peak_near) = group_num;
            peak_pool(peak_near) = 0;
            
            % Set newest layer member as the next peak to use in search
            peak_n = peak_near;
            
        else
            % If no neighbors are within tolerance, end search to left of
            % layer
            search_L = false;
        end
    end
    
%     layer_i = find(peak_group==group_num);
%     peak_pool(layer_i) = 0;
    

    % Find subscript indices of all members of ith layer
    [row, col] = ind2sub(size(peaks_raw), layer_i);
    
    % Smooth layer using a moving mean of row positions
    row = round(movmean(row, 10));
    
    % For loop to search for additional peaks near current layer
    for j = 1:length(layer_i)
        
        % Find 5 nearest layer members on either side of jth member
        mag_idx = sub2ind(size(peaks_raw), ...
            row(max([1 j-5]):min([length(row) j+5])), ...
            col(max([1 j-5]):min([length(col) j+5])));
        
        % Find local mean peak magnitude and width from the nearby layer
        % members
        mag_j = mean(peaks_raw(mag_idx));
        width_j = mean(peak_width(mag_idx));
        
        % Define local search window as 0.50 m above and below jth member
        % and 100 m to left and right of jth member
        row_idx = [max([row(j)-round(0.50/core_res) 1]) ...
            min([row(j)+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col(j)-round(100/horz_res)]) ...
            min([size(peaks_raw, 2) col(j)+round(100/horz_res)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        % Calculate distance between jth layer member and peaks within the
        % local search window
        dist_j = sqrt((row(j) - (row_idx(1)+row_local-1)).^2 + ...
            2*(col(j) - (col_idx(1)+col_local-1)).^2 + (mag_j - mag_local).^2);
        
        % Set distance threshold based on peak width and error bin size
        threshold = 1*width_j + err_bin;
        
        % Assign all peaks within tolerance to the current layer group and
        % remove those peaks from the pool of available peaks to search
        tol_idx = dist_j <= threshold;
        group_idx = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(tol_idx)-1, col_idx(1)+col_local(tol_idx)-1);
        peak_group(group_idx) = group_num;
        peak_pool(group_idx) = 0;
        
    end
    
    % Find all members of the ith layer group and output to preallocated
    % array
    layer_i = matrix_idx(peak_group==group_num);
    layers{i} = layer_i;
    
    % Initialize values for next iteration of while loop
    group_num = group_num + 1;
    i = i + 1;
    
    % End search for new layers when remaining ungrouped peaks are smaller 
    % than a given threshold
    if max(peak_pool(:)) <= minProm
        search_new = false;
    end
end











% % While loop that iterates on each accumulation layer
% while search_new == true
%     
%     % Find the maximum peak remaining in pool and assign to group
%     [~, peak_max] = max(peak_pool(:));
%     peak_group(peak_max) = group_num;
%     
%     % Initialize values for group search within layer to right of max peak
%     peak_n = peak_max;
%     search_R = true;
%     
%     % Search to right of most recent group member for new group members
%     while search_R == true
%         
%         % Find all current layer group members
%         layer_i = matrix_idx(peak_group==group_num);
%         
%         % Get nearest 5 group members
%         group_idx = layer_i(max([1 length(layer_i)-4]):length(layer_i));
%         
%         % Find subscripts of farthest right layer group member
%         [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
%         
%         % Define estimated row position for projected nearest neighbor as 
%         % the weighted average of nearest 5 layer members and the farthest
%         % right layer member (50/50 weighting)
%         [row_i, ~] = ind2sub(size(peaks_raw), group_idx);
%         row_n = round((1/2)*(row_n + mean(row_i)));
%         
%         % Define estimated peak magnitude and peak width for the projected  
%         % nearest neighbor as the mean magnitude and width of the nearest 
%         % 5 group members
%         mag_n = mean(peaks_raw(group_idx));
%         width_n = mean(peak_width(group_idx));
%         
%         % Define local  search window as 250 m to right of last known group
%         % member, and 0.50 m above and below estimated row position for the
%         % projected nearest neighbor
%         row_idx = [max([row_n-round(0.50/core_res) 1]) ...
%             min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
%         col_idx = [col_n+1 min([col_n+round(250/horz_res) size(peaks_raw, 2)])];
%         
%         % Get matrix of peaks from avaiable pool within the local search
%         % area, along with their indices, magnitudes, and subscript indices
%         peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
%         local_idx = find(peak_local);
%         mag_local = peak_local(local_idx);
%         [row_local, col_local] = ind2sub(size(peak_local), local_idx);
%         
%         % Calculate distance between projected layer position and peaks
%         % within the local search window
%         dist_n = sqrt((row_n - (row_idx(1)+row_local-1)).^2 + ...
%             2*(col_n - (col_idx(1)+col_local-1)).^2 + (mag_n - mag_local).^2);
%         
%         % Select the nearest neighbor to the estimated layer position
%         [min_dist, dist_idx] = min(dist_n);
%         
%         % Set distance threshold based on peak width and error bin size
%         threshold = 1*width_n + err_bin;
%         
%         % Determine if nearest neighbor is within distance tolerance
%         if min_dist <= threshold
%             
%             % Find index of nearest neighbor (when within tolerance)
%             peak_near = sub2ind(size(peaks_raw), ...
%                 row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
%             
%             % Assign nearest neighbor to current layer group and remove
%             % from the pool of available peaks to search
%             peak_group(peak_near) = group_num;
%             peak_pool(peak_near) = 0;
%             
%             % Set newest layer member as the next peak to use in search
%             peak_n = peak_near;
%             
%         else
%             % If no neighbors are within tolerance, end search to right of
%             % layer
%             search_R = false;
%         end
%     end
%     
%     % Initialize values for group search within layer to left of max peak
%     peak_n = peak_max;
%     search_L = true;
%     
%     % Search to left of most recent group member for new group members
%     while search_L == true
%         
%         % Find all current layer group members
%         layer_i = matrix_idx(peak_group==group_num);
%         
%         % Get nearest 5 group members
%         group_idx = layer_i(1:min([5 length(layer_i)]));
%         
%         % Find subscripts of farthest left layer member
%         [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
%         
%         % Define estimated row position for projected nearest neighbor as 
%         % the weighted average of nearest 5 layer members and the farthest
%         % left layer member (50/50 weighting)
%         [row_i, ~] = ind2sub(size(peaks_raw), group_idx);
%         row_n = round((1/2)*(row_n + mean(row_i)));
%         
%         % Define estimated peak magnitude and peak width for the projected  
%         % nearest neighbor as the mean magnitude and width of the nearest 
%         % 5 group members
%         mag_n = mean(peaks_raw(group_idx));
%         width_n = mean(peak_width(group_idx));
%         
%         % Define local  search window as 250 m to left of last known group
%         % member, and 0.50 m above and below estimated row position for the
%         % projected nearest neighbor
%         row_idx = [max([row_n-round(0.50/core_res) 1]) ...
%             min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
%         col_idx = [max([1 col_n-round(250/horz_res)-1]) col_n-1];
%         
%         % Get matrix of peaks from avaiable pool within the local search
%         % area, along with their indices, magnitudes, and subscript indices
%         peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
%         local_idx = find(peak_local);
%         mag_local = peak_local(local_idx);
%         [row_local, col_local] = ind2sub(size(peak_local), local_idx);
%         
%         % Calculate distance between projected layer position and peaks
%         % within the local search window
%         dist_n = sqrt((row_n - (row_idx(1)+row_local-1)).^2 + ...
%             2*(col_n - (col_idx(1)+col_local-1)).^2 + (mag_n - mag_local).^2);
%         
%         % Select the nearest neighbor to the estimated layer position
%         [min_dist, dist_idx] = min(dist_n);
%         
%         % Set distance threshold based on peak width and error bin size
%         threshold = 1*width_n + err_bin;
%         
%         % Determine if nearest neighbor is within distance tolerance
%         if min_dist <= threshold
%             
%             % Find index of nearest neighbor (when within tolerance)
%             peak_near = sub2ind(size(peaks_raw), ...
%                 row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
%             
%             % Assign nearest neighbor to current layer group and remove
%             % from the pool of available peaks to search
%             peak_group(peak_near) = group_num;
%             peak_pool(peak_near) = 0;
%             
%             % Set newest layer member as the next peak to use in search
%             peak_n = peak_near;
%             
%         else
%             % If no neighbors are within tolerance, end search to left of
%             % layer
%             search_L = false;
%         end
%     end
%     
% %     layer_i = find(peak_group==group_num);
% %     peak_pool(layer_i) = 0;
%     
% 
%     % Find subscript indices of all members of ith layer
%     [row, col] = ind2sub(size(peaks_raw), layer_i);
%     
%     % Smooth layer using a moving mean of row positions
%     row = round(movmean(row, 10));
%     
%     % For loop to search for additional peaks near current layer
%     for j = 1:length(layer_i)
%         
%         % Find 5 nearest layer members on either side of jth member
%         mag_idx = sub2ind(size(peaks_raw), ...
%             row(max([1 j-5]):min([length(row) j+5])), ...
%             col(max([1 j-5]):min([length(col) j+5])));
%         
%         % Find local mean peak magnitude and width from the nearby layer
%         % members
%         mag_j = mean(peaks_raw(mag_idx));
%         width_j = mean(peak_width(mag_idx));
%         
%         % Define local search window as 0.50 m above and below jth member
%         % and 100 m to left and right of jth member
%         row_idx = [max([row(j)-round(0.50/core_res) 1]) ...
%             min([row(j)+round(0.50/core_res) size(peaks_raw, 1)])];
%         col_idx = [max([1 col(j)-round(100/horz_res)]) ...
%             min([size(peaks_raw, 2) col(j)+round(100/horz_res)])];
%         
%         % Get matrix of peaks from avaiable pool within the local search
%         % area, along with their indices, magnitudes, and subscript indices
%         peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
%         local_idx = find(peak_local);
%         mag_local = peak_local(local_idx);
%         [row_local, col_local] = ind2sub(size(peak_local), local_idx);
%         
%         % Calculate distance between jth layer member and peaks within the
%         % local search window
%         dist_j = sqrt((row(j) - (row_idx(1)+row_local-1)).^2 + ...
%             2*(col(j) - (col_idx(1)+col_local-1)).^2 + (mag_j - mag_local).^2);
%         
%         % Set distance threshold based on peak width and error bin size
%         threshold = 1*width_j + err_bin;
%         
%         % Assign all peaks within tolerance to the current layer group and
%         % remove those peaks from the pool of available peaks to search
%         tol_idx = dist_j <= threshold;
%         group_idx = sub2ind(size(peaks_raw), ...
%             row_idx(1)+row_local(tol_idx)-1, col_idx(1)+col_local(tol_idx)-1);
%         peak_group(group_idx) = group_num;
%         peak_pool(group_idx) = 0;
%         
%     end
%     
%     % Find all members of the ith layer group and output to preallocated
%     % array
%     layer_i = matrix_idx(peak_group==group_num);
%     layers{i} = layer_i;
%     
%     % Initialize values for next iteration of while loop
%     group_num = group_num + 1;
%     i = i + 1;
%     
%     % End search for new layers when remaining ungrouped peaks are smaller 
%     % than a given threshold
%     if max(peak_pool(:)) <= minProm
%         search_new = false;
%     end
% end

% Remove empty cells and layers shorter than 100 m
layers = layers(cellfun(@(x) length(x) > 4, layers));

%%

% Preallocate arrays for the matrix indices of members of each layer
layers_idx = cell(1,length(layers));
% layers_test = cell(1,length(layers));
peaks = zeros(size(peaks_raw));

% For loop to coerce layers to have one row position for each trace
for i = 1:length(layers_idx)
    
    % Find matrix indices of all members of ith layer
    layer_i = layers{i};
    
    % Find row and col indices of members of ith layer
    [row, col] = ind2sub(size(radar.data_smooth), layer_i);
    
%     [~,sort_idx] = sort(col);
%     row = row(sort_idx);
%     col = col(sort_idx);
    
    % If multiple rows exist for the same column, take the
    % squared prominence-weighted mean of the rows
    if length(col) > length(unique(col))
        
        % Create matrix of peaks with only layer_i members
        layer_mat = zeros(size(radar.data_smooth));
        layer_mat(layer_i) = peaks_raw(layer_i);
        
        % Find all traces with multiple peaks within layer_i
        multi_idx = sum(logical(layer_mat))>1;
        col_nums = 1:size(layer_mat, 2);
        for k = col_nums(multi_idx)
            k_idx = find(col==k);
            k_peaks = layer_mat(row(k_idx),k);
            
            % Create squared prominence weighting matrix
            k_sum = sum(k_peaks.^2);
            k_row = round(sum((k_peaks.^2/k_sum).*row(k_idx)));
            k_peak = sum((k_peaks.^2/k_sum).*k_peaks);
%             k_sum = sum(k_peaks);
%             k_row = round(sum((k_peaks/k_sum).*row(k_idx)));
%             k_peak = sum((k_peaks/k_sum).*k_peaks);

            % Assign squared prominence weighted mean position to layer
            % trace
            k_col = zeros(size(layer_mat, 1), 1);
            k_col(k_row) = k_peak;
            layer_mat(:,k) = k_col;
        end
        peaks_mat = find(layer_mat);
        [row, col] = ind2sub(size(radar.data_smooth), peaks_mat);
        peaks(peaks_mat) = layer_mat(peaks_mat);
    
    else
        % Assign weighted peak position and value to preallocated peaks
        % matrix
        peaks(layer_i) = peaks_raw(layer_i);
    end

    % Smooth layer i using a moving average of row indices and assign 
    % smoothed layer to preallocated cell array
    row_mean = round(movmean(row, round(100/horz_res)));
    c_interp = min(col):max(col);
    r_interp = round(interp1(col, row_mean, c_interp));
    layers_idx{i} = sub2ind(size(radar.data_smooth), r_interp, c_interp);
%     layers_idx{i} = sub2ind(size(radar.data_smooth), row_mean, col);
%     layers_idx{i} = sub2ind(size(radar.data_smooth), row, col);
    
end

% Calculate continuous layer distances for each layer (accounting for 
% lateral size of stacked radar trace bins)
% layers_dist = cellfun(@(x) numel(x)*horz_res, layers_idx);
layers_dist = cellfun(@(x) numel(x)*horz_res, layers_idx);

% Map layer prominence-distance values to the location within the radar
% matrix of the ith layer
layer_peaks = zeros(size(peaks));
for i = 1:length(layers_idx)
    layer_peaks(layers_idx{i}) = sum(peaks(layers_idx{i})).*layers_dist(i);
end


% Output layer arrays to radar structure
radar.peaks = peaks;
radar.layers = layers_idx;
radar.layer_vals = layer_peaks;

%% Assign layer likelihood scores and estimate age-depth scales

ages = zeros([size(radar.data_smooth) Ndraw]);
radar.likelihood = zeros(size(radar.data_smooth));
for i = 1:size(layer_peaks, 2)
    
%     P_50 = 2*1000;
%     P_50 = 1000*mean(std(radar.data_smooth));
%     P_50 = 1000*(quantile(radar.data_smooth(:,i), 0.95) - ...
%         quantile(radar.data_smooth(:,i), 0.05));
    P_50 = 1*500*median(Proms{i});
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
