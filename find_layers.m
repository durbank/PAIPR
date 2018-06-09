function [peak_group, layers] = find_layers(peaks_raw, peak_width, err_bin, core_res, horz_res)

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
    layer_i = matrix_idx(peak_group==group_num);
    search_R = true;
    
    while search_R == true && length(layer_i) < 4
        
        % Get all group members
        group_idx = layer_i;
        
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
        [~, dist_idx] = min(dist_n);
        
        
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
        if (col_idx(1)+col_local(dist_idx)-1) >= size(peaks_raw, 2)
            search_R = false;
        end
        
        % Update list of current group members
        layer_i = matrix_idx(peak_group==group_num);
    end
    
    
    peak_n = peak_max;
    
    while length(layer_i) < 4
        
        % Get all group members
        group_idx = layer_i;
        
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
        
        % Define local  search window as 250 m to left of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_n-round(0.50/core_res) 1]) ...
            min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col_n-round(250/horz_res)-1]) max([col_n-1 1])];
        
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
        [~, dist_idx] = min(dist_n);
        
        % Find index of nearest neighbor (when within tolerance)
        peak_near = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
        
        % Assign nearest neighbor to current layer group and remove
        % from the pool of available peaks to search
        peak_group(peak_near) = group_num;
        peak_pool(peak_near) = 0;
        
        % Set newest layer member as the next peak to use in search
        peak_n = peak_near;
        
        %         % If layer has reached the right edge of radargram, end search
        %         % to right of layer
        %         if (col_idx(1)+col_local(min_idx)-1) >= size(peaks_raw, 2)
        %             search_L = false;
        %         end
        
        % Find all current layer group members
        layer_i = matrix_idx(peak_group==group_num);
    end
    
    
    % Initialize values for group search within layer to right of max peak
    peak_n = layer_i(end);
    search_R = true;
    
    while search_R == true
        
        % Get row, col, and magnitude of current members of layer group
        [row_i, col_i] = ind2sub(size(peaks_raw), layer_i);
        mag_i = peaks_raw(layer_i);
        
        % Get nearest 10 group members
        group_idx = layer_i(end-min([10 length(layer_i)])+1:end);
        
        % Find subscripts of farthest right layer group member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        
        % Define estimated row position for projected nearest neighbor as
        % the weighted average of nearest 5 layer members and the farthest
        % right layer member (50/50 weighting)
        [row_temp, ~] = ind2sub(size(peaks_raw), group_idx);
        row_n = round((1/2)*(row_n + mean(row_temp)));
        col_n = col_n + 1;
        
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
        
        % MAY NEED TO CHANGE TO LOCAL GROUP MEMBERS, NOT ALL GROUP MEMBERS
        data_i = [row_i col_i mag_i];
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1) mag_local];
        
        data_EX = repmat([row_n col_n mag_n], length(local_idx), 1);
        
        
        sigma = cov(data_i);
        
        dist_n = diag((data_local - data_EX)*inv(sigma)*(data_local - data_EX)');
        
        % Select the nearest neighbor to the estimated layer position
        [min_dist, dist_idx] = min(dist_n);
        
        % Set distance threshold based on peak width and error bin size
        threshold = 100;
        
        % ADD ALL LOCAL PEAKS WITHIN THRESHOLD TO GROUP, NOT JUST NEAREST
        
        
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
            
            % Update list of current group members
            layer_i = matrix_idx(peak_group==group_num);
            % layer_i = sort([layer_i; peak_near]);
            
            % If layer has reached the right edge of radargram, end search
            % to right of layer
            if (col_idx(1)+col_local(dist_idx)-1) >= size(peaks_raw, 2)
                search_R = false;
            end
            
        else
            % If no neighbors are within tolerance, end search to right of
            % layer
            search_R = false;
            
            
        end
        
    end
    
    
    
    
    
    
    
    
    
end
end