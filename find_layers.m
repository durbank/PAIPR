function [peak_group, layers] = find_layers(peaks_raw, peak_width, core_res, horz_res)

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
    
%         A = peak_pool;
%         B = ones(6)/6^2;
%         C = conv2(A,B,'same');
%         [~,C_max] = max(C(:));
%         [r_C, c_C] = ind2sub(size(C), C_max);
    
    % Find the maximum peak remaining in pool and assign to group
    
    [~, peak_max] = max(peak_pool(:));
    peak_group(peak_max) = group_num;
    peak_pool(peak_max) = 0;
    
    % Initialize values for group search within layer to right of max peak
    peak_n = peak_max;
    layer_i = matrix_idx(peak_group==group_num);
    search_R = true;
    search_L = true;
    
    while search_R == true && length(layer_i) < 5
        
        % Get all group members
        group_idx = layer_i;
        
        % Find subscripts of farthest right layer group member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        
        % Define estimated row position for projected nearest neighbor as
        % the weighted average of nearest 5 layer members and the farthest
        % right layer member (50/50 weighting)
        [row_temp, ~] = ind2sub(size(peaks_raw), group_idx);
        row_pred = round(median(row_temp));
        col_pred = col_n + 1;
        
        % Define estimated peak magnitude and peak width for the projected
        % nearest neighbor as the mean magnitude and width of the nearest
        % 5 group members
        mag_n = median(peaks_raw(group_idx));
        width_n = median(peak_width(group_idx));
        
        % Define local  search window as 250 m to right of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_pred-round(0.50/core_res) 1]) ...
            min([row_pred+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [min([col_pred size(peaks_raw, 2)]) ...
            min([col_pred+round(250/horz_res) size(peaks_raw, 2)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1) mag_local];
        
        % Calculate distance between projected layer position and peaks
        % within the local search window
        dist_n = sqrt(((data_local(:,1)-row_pred)/(0.5*width_n)).^2 + ...
            ((data_local(:,2)-col_pred)/(0.5*(col_idx(2)-col_idx(1)))).^2 + ...
            ((mag_local-mag_n)/1).^2);
        %         dist_n = sqrt((row_pred - (row_idx(1)+row_local-1)).^2 + ...
        %             2*(col_pred - (col_idx(1)+col_local-1)).^2 + (mag_n - mag_local).^2);
        
        % Select the nearest neighbor to the estimated layer position
        %         [~, dist_idx] = min(dist_n);
        %         threshold = median(peak_width(group_idx)) + err_bin;
        threshold = 3;
        dist_idx = dist_n <= threshold;
        
        
        
        % Find index of nearest neighbor (when within tolerance)
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
        
        % Assign nearest neighbor to current layer group and remove
        % from the pool of available peaks to search
        peak_group(peak_neighbor) = group_num;
        peak_pool(peak_neighbor) = 0;
        
        % Update list of current group members
        layer_i = matrix_idx(peak_group==group_num);
        
        % Set newest layer member as the next peak to use in search
        peak_n = layer_i(end);
        
        if isempty(peak_neighbor)
            search_R = false;
        end
        
        % If layer has reached the right edge of radargram, end search
        % to right of layer
        if (max(col_idx(1)+col_local(dist_idx)-1)) >= size(peaks_raw, 2)
            search_R = false;
        end
        
        
    end
    
    
    peak_n = peak_max;
    
    while length(layer_i) < 5 && search_L == true
        
        % Get all group members
        group_idx = layer_i;
        
        % Find subscripts of farthest left layer group member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        
        % Define estimated row position for projected nearest neighbor as
        % the weighted average of nearest 5 layer members and the farthest
        % right layer member (50/50 weighting)
        [row_temp, ~] = ind2sub(size(peaks_raw), group_idx);
        row_pred = round(median(row_temp));
        %         row_n = round((1/2)*(row_n + median(row_temp)));
        col_pred = col_n - 1;
        
        % Define estimated peak magnitude and peak width for the projected
        % nearest neighbor as the mean magnitude and width of the nearest
        % 5 group members
        mag_n = median(peaks_raw(group_idx));
        width_n = median(peak_width(group_idx));
        
        % Define local  search window as 250 m to left of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_pred-round(0.50/core_res) 1]) ...
            min([row_pred+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col_pred-round(150/horz_res)]) max([col_pred 1])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1) mag_local];
        
        % Calculate distance between projected layer position and peaks
        % within the local search window
        dist_n = sqrt(((data_local(:,1)-row_pred)/(0.5*width_n)).^2 + ...
            ((data_local(:,2)-col_pred)/(0.5*(col_idx(2)-col_idx(1)))).^2 + ...
            ((mag_local-mag_n)/1).^2);
        %         dist_n = sqrt((row_n - (row_idx(1)+row_local-1)).^2 + ...
        %             2*(col_n - (col_idx(1)+col_local-1)).^2 + (mag_n - mag_local).^2);
        
        % Select the nearest neighbor to the estimated layer position
        [~, dist_idx] = min(dist_n);
        
        % Find index of nearest neighbor (when within tolerance)
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
        
        if isempty(peak_neighbor)
            search_L = false;
            
        else
            % Assign nearest neighbor to current layer group and remove
            % from the pool of available peaks to search
            peak_group(peak_neighbor) = group_num;
            peak_pool(peak_neighbor) = 0;
            
            % Set newest layer member as the next peak to use in search
            peak_n = peak_neighbor;
            
            % If layer has reached the right edge of radargram, end search
            % to right of layer
            if (col_idx(1)+col_local(dist_idx)-1) <= 1
                search_L = false;
            end
            
            % Find all current layer group members
            layer_i = matrix_idx(peak_group==group_num);
        end
        
    end
    
    
    % Initialize values for group search within layer to right of max peak
    %     peak_n = layer_i(end);
    if search_L == false
        search_R = false;
    else
        search_R = true;
    end
    
    while search_R == true
        
        % Get nearest 10 group members
        group_idx = layer_i(end-min([19 length(layer_i)])+1:end);
        
        % Get row, col, and magnitude of current members of layer group
        [row_i, col_i] = ind2sub(size(peaks_raw), group_idx);
        mag_i = peaks_raw(group_idx);
        
        % Find subscripts of farthest right layer group member
        %         [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        row_n = row_i(end);
        col_n = col_i(end);
        
        % Define estimated row position for projected nearest neighbor as
        % the weighted average of nearest 10 layer members and the farthest
        % right layer member (50/50 weighting)
        %         [row_temp, ~] = ind2sub(size(peaks_raw), group_idx);
        %%USE FULL WEIGHTED ROW MEAN IN THE FUTURE
        %         row_pred = round(median(row_i(end-4:end)));
        col_pred = col_n + 1;
        %         row_bar = movmedian(row_i, 4);
        %         p = polyfit(col_i, row_i, 1);
        %         row_pred = round(polyval(p, col_pred));
        try
            pp = csaps(col_i, movmedian(row_i, 4), 0.05);
            row_pred = round(fnval(pp, col_pred));
        catch
            row_pred = round(median(row_i(end-4:end)));
        end
        
        %         row_pred = round(median(row_i(end-4:end)));
        
        % Define estimated peak magnitude and peak width for the projected
        % nearest neighbor as the mean magnitude and width of the nearest
        % 5 group members
        mag_n = median(peaks_raw(group_idx));
        mag_var = var(peaks_raw(group_idx));
        width_n = median(peak_width(group_idx));
        
        % Define local  search window as 250 m to right of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_pred-round(0.50/core_res) 1]) ...
            min([row_pred+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [col_pred min([col_pred+round(150/horz_res) size(peaks_raw, 2)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        % MAY NEED TO CHANGE TO LOCAL GROUP MEMBERS, NOT ALL GROUP MEMBERS
        data_i = [row_i col_i mag_i];
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1) mag_local];
        data_EX = repmat([row_pred col_pred mag_n], length(local_idx), 1);
        sigma = cov(data_i);
        
        dist_n = sqrt(((data_local(:,1)-row_pred)/(0.5*width_n)).^2 + ...
            ((data_local(:,2)-col_pred)/(0.5*(col_idx(2)-col_idx(1)))).^2 + ...
            ((mag_local-mag_n)/1).^2);
        %         dist_n = sqrt(diag((data_local - data_EX)*inv(sigma)*...
        %             (data_local - data_EX)'));
        %         dist_n = pdist2(data_local, data_EX(1,:), 'mahalanobis', sigma);
        
        % Select the nearest neighbor to the estimated layer position
        %         [min_dist, dist_idx] = min(dist_n);
        
        % Set distance threshold (based on p <= 0.01 for n-1 df on
        % Chi-squared distribution)
        %         threshold = 5.991;
        threshold = 3;
        dist_idx = dist_n <= threshold;
        
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
        
        peak_group(peak_neighbor) = group_num;
        peak_pool(peak_neighbor) = 0;
        
        % Update list of current group members
        layer_i = matrix_idx(peak_group==group_num);
        
        % Set newest layer member as the next peak to use in search
        %         peak_n = layer_i(end);
        
        % If layer has reached the right edge of radargram, end search
        % to right of layer
        if (max(col_idx(1)+col_local(dist_idx)-1)) >= size(peaks_raw, 2)
            search_R = false;
        end
        
        if isempty(peak_neighbor)
            search_R = false;
            
            % ADD ALL LOCAL PEAKS WITHIN THRESHOLD TO GROUP, NOT JUST NEAREST
        end
        
    end
    
    % Initialize values for group search within layer to right of max peak
    %     peak_n = layer_i(1);
    %     search_L = true;
    
    while search_L == true
        
        % Get nearest 10 group members
        group_idx = layer_i(1:min([20 length(layer_i)]));
        
        % Get row, col, and magnitude of current members of layer group
        [row_i, col_i] = ind2sub(size(peaks_raw), group_idx);
        mag_i = peaks_raw(group_idx);
        
        % Find subscripts of farthest right layer group member
        %         [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        row_n = row_i(1);
        col_n = col_i(1);
        
        % Define estimated row position for projected nearest neighbor as
        % the weighted average of nearest 10 layer members and the farthest
        % right layer member (50/50 weighting)
        %         [row_temp, ~] = ind2sub(size(peaks_raw), group_idx);
        %%USE FULL WEIGHTED ROW MEAN IN THE FUTURE
        col_pred = col_n - 1;
        
        try
            pp = csaps(col_i, movmedian(row_i, 4), 0.05);
            row_pred = round(fnval(pp, col_pred));
        catch
            row_pred = round(median(row_i(1:5)));
        end
        %         pp = csaps(col_i, row_i, 0.05);
        %         row_pred = round(fnval(pp, col_pred));
        %         p = polyfit(col_i, row_i, 1);
        %         row_pred = round(polyval(p, col_pred));
        %         row_pred = round(median(row_i(1:5)));
        
        % Define estimated peak magnitude and peak width for the projected
        % nearest neighbor as the mean magnitude and width of the nearest
        % 20 group members
        mag_n = median(peaks_raw(group_idx));
        mag_var = var(peaks_raw(group_idx));
        width_n = median(peak_width(group_idx));
        
        % Define local  search window as 250 m to right of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_pred-round(0.50/core_res) 1]) ...
            min([row_pred+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col_pred-round(150/horz_res)]) col_pred];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        
        % MAY NEED TO CHANGE TO LOCAL GROUP MEMBERS, NOT ALL GROUP MEMBERS
        data_i = [row_i col_i mag_i];
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1) mag_local];
        data_EX = repmat([row_pred col_pred mag_n], length(local_idx), 1);
        sigma = cov(data_i);
        
        dist_n = sqrt(((data_local(:,1)-row_pred)/(0.5*width_n)).^2 + ...
            ((data_local(:,2)-col_pred)/(0.5*(col_idx(2)-col_idx(1)))).^2 + ...
            ((mag_local-mag_n)/1).^2);
        %         dist_n = sqrt(diag((data_local - data_EX)*inv(sigma)*...
        %             (data_local - data_EX)'));
        %         dist_n = pdist2(data_local, data_EX(1,:), 'mahalanobis', sigma);
        
        % Select the nearest neighbor to the estimated layer position
        %         [min_dist, dist_idx] = min(dist_n);
        
        % Set distance threshold (based on p <= 0.01 for n-1 df on
        % Chi-squared distribution)
        %         threshold = 5.991;
        threshold = 3;
        dist_idx = dist_n <= threshold;
        
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
        
        peak_group(peak_neighbor) = group_num;
        peak_pool(peak_neighbor) = 0;
        
        % Update list of current group members
        layer_i = matrix_idx(peak_group==group_num);
        
        % Set newest layer member as the next peak to use in search
        %         peak_n = layer_i(1);
        
        % If layer has reached the right edge of radargram, end search
        % to right of layer
        if (min(col_idx(1)+col_local(dist_idx)-1)) <= 1
            search_L = false;
        end
        
        if isempty(peak_neighbor)
            search_L = false;
        end
        
    end
    
    
    
    
    
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
    
    
    
    
    % Find all members of the ith layer group and output to preallocated
    % array
    layers{i} = layer_i;
    
    % Initialize values for next iteration of while loop
    group_num = group_num + 1;
    i = i + 1;
    
    % End search for new layers when remaining ungrouped peaks are smaller
    % than a given threshold
    if max(peak_pool(:)) <= 0.50
        search_new = false;
    end
    
end

% Remove empty cells and layers shorter than 100 m
layers = layers(cellfun(@(x) length(x) > 4, layers));

end