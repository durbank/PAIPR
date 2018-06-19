% Function that builds continuous layers by adding only nearest neighbor
% (rather than all that meet threshold) based on the position of last
% member added (rather than predicted position). The method uses a modified
% weighted Euclidean distance function (row value is scaled by the peak 
% half-width) to find the nearest neighbor

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
    
    % Convolve remaining peaks with a 10x10 matrix and find brightest 
    % convolved pixel (2D smoothing of matrix in order to find brightest 
    % local region)
    A = peak_pool;
    B = ones(10)/10^2;
    C = conv2(A,B,'same');
    [~,C_max] = max(C(:));
    [r_C, c_C] = ind2sub(size(C), C_max);
    
    % Find local area around max convolved peak
    r_idx = [max([1 r_C-5]) min([size(peaks_raw,1) r_C+5])];
    c_idx = [max([1 c_C-5]) min([size(peaks_raw,2) c_C+5])];
    
    % Find index of actual brightest peak within local search area
    peak_local = peak_pool(r_idx(1):r_idx(2),c_idx(1):c_idx(2));
    [~, local_max] = max(peak_local(:));
    [r_max, c_max] = ind2sub(size(peak_local), local_max);
    peak_max = sub2ind(size(peak_pool), r_idx(1)+r_max-1, c_idx(1)+c_max-1);
    

    % Assign max peak to current group and remove from future peak searches
%     [~, peak_max] = max(peak_pool(:));
    peak_group(peak_max) = group_num;
    peak_pool(peak_max) = 0;
    
    % Initialize values for group search within layer to right of max peak
    peak_n = peak_max;
    layer_i = matrix_idx(peak_group==group_num);
    search_R = true;
    
    
    while search_R == true
        
        % Get nearest 20 group members
        group_idx = layer_i(end-min([19 length(layer_i)])+1:end);
        
        % Get median magnitude and median peak width of 20 nearest group
        % members
        mag_i = median(peaks_raw(group_idx));
        width_i = median(peak_width(group_idx));
        
%         if length(group_idx) <= 4
%             mag_var = 1;
%             
%         else
%             mag_var = var(peaks_raw(group_idx));
%         end
        
        % Find subscripts of farthest right layer group member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        col_n = col_n + 1;
        
        
        % Define local  search window as 150 m to right of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_n-round(0.50/core_res) 1]) ...
            min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [col_n min([col_n+round(150/horz_res) size(peaks_raw, 2)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1)];
        
        % Calculate the weighted Euclidean distance between last known
        % layer group member and peaks within the local window, where peak
        % magnitude is weighted 0.5 and the row distance is scaled by the
        % layer's half-width
        dist_n = sqrt(1*(((data_local(:,1)-row_n)/(0.5*width_i))).^2 + ...
            1*(data_local(:,2)-col_n).^2 + 0.5*(mag_local-mag_i).^2);
        
        % Set distance threshold for group membership and find the nearest
        % neighbor within local search
        threshold = 4;
        [min_dist, dist_idx] = min(dist_n);
        
        % Determine if nearest neighbor is within tolerance of layer group
        if min_dist <= threshold
            
            % Get matrix index of nearest neighbor
            peak_neighbor = sub2ind(size(peaks_raw), ...
                row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
            
            % Assign nearest neighbor to layer group, and remove peak from
            % future searches
            peak_group(peak_neighbor) = group_num;
            peak_pool(peak_neighbor) = 0;
            
            % Update list of current group members
            peak_n = peak_neighbor;
            layer_i = matrix_idx(peak_group==group_num);
            
            % Check if newest group member is in the final data column to
            % the right
            [~, col_check] = ind2sub(size(peaks_raw), peak_n);
            
            % If layer has reached the right edge of radargram, end search
            % to right of layer
            if col_check >= size(peaks_raw, 2)
                search_R = false;
            end
            
        else
            % If nearest right neighbor is outside group tolerance, end
            % search to right of layer
            search_R = false;
        end
        
    end
    
    
    % Initialize values for search to left of max peak
    peak_n = peak_max;
    search_L = true;
    
    while search_L == true
        
        % Get nearest 20 group members
        group_idx = layer_i(1:min([20 length(layer_i)]));
        
        % Get median magnitude and median peak width of 20 nearest group
        % members
        mag_i = median(peaks_raw(group_idx));
        width_i = median(peak_width(group_idx));
        
%         if length(group_idx) <= 4
%             mag_var = 1;
%             
%         else
%             mag_var = var(peaks_raw(group_idx));
%         end
        
        
        
        % Find subscripts of farthest left layer group member
        [row_n, col_n] = ind2sub(size(peaks_raw), peak_n);
        col_n = col_n - 1;
        
        
        % Define local  search window as 150 m to left of last known group
        % member, and 0.50 m above and below most recent group member
        row_idx = [max([row_n-round(0.50/core_res) 1]) ...
            min([row_n+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col_n-round(150/horz_res)]) col_n];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1) mag_local];
        
        % Calculate the weighted Euclidean distance between last known
        % layer group member and peaks within the local window, where peak
        % magnitude is weighted 0.5 and the row distance is scaled by the
        % layer's half-width
        dist_n = sqrt(1*(((data_local(:,1)-row_n)/(0.5*width_i))).^2 + ...
            1*(data_local(:,2)-col_n).^2 + 0.5*(mag_local-mag_i).^2);
        
        % Set distance threshold for group membership and find the nearest
        % neighbor within local search
        threshold = 4;
        [min_dist, dist_idx] = min(dist_n);
        
        % Determine if nearest neighbor is within tolerance of layer group
        if min_dist <= threshold
            
            % Get matrix index of nearest neighbor
            peak_neighbor = sub2ind(size(peaks_raw), ...
                row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
            
            % Assign nearest neighbor to layer group, and remove peak from
            % future searches
            peak_group(peak_neighbor) = group_num;
            peak_pool(peak_neighbor) = 0;
            
            % Update list of current group members
            peak_n = peak_neighbor;
            layer_i = matrix_idx(peak_group==group_num);
            
            % Check if newest group member is in the final data column to
            % the right
            [~, col_check] = ind2sub(size(peaks_raw), peak_n);
            
            % If layer has reached the left edge of radargram, end search
            % to left of layer
            if col_check <= 1
                search_L = false;
            end
            
        else
            % If nearest left neighbor is outside group tolerance, end
            % search to left of layer
            search_L = false;
        end
    end
    

    
    % Find subscript indices of all members of ith layer
    [row, col] = ind2sub(size(peaks_raw), layer_i);
    
    % Smooth layer using locally weighted regression
%     row_mean = round(movmean(row, 5));
    row_mean = round(smooth(col, row, 'rlowess'));
    
    % For loop to search for additional peaks near current layer
    for j = 1:length(layer_i)
        
        % Find 10 nearest layer members on either side of jth member
        group_idx = sub2ind(size(peaks_raw), ...
            row(max([1 j-10]):min([length(row) j+10])), ...
            col(max([1 j-10]):min([length(col) j+10])));
        
        % Find local mean peak magnitude and width from the nearby layer
        % members
        mag_j = median(peaks_raw(group_idx));
        width_j = median(peak_width(group_idx));
        
        % Define local search window as 0.50 m above and below jth member
        % and 150 m to left and right of jth member
        row_idx = [max([row_mean(j)-round(0.50/core_res) 1]) ...
            min([row_mean(j)+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col(j)-round(150/horz_res)]) ...
            min([size(peaks_raw, 2) col(j)+round(150/horz_res)])];
        
        % Get matrix of peaks from available pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1)];
        
        % Calculate distance between jth layer member and peaks within the
        % local search window
        dist_j = sqrt(((data_local(:,1)-row_mean(j))/(0.5*width_j)).^2 + ...
            (data_local(:,2)-col(j)).^2 + 0.5*(mag_local - mag_j).^2);
        
        % Set distance threshold
        threshold = 3;
        
        % Assign all peaks within tolerance to the current layer group and
        % remove those peaks from the pool of available peaks to search
        tol_idx = dist_j <= threshold;
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(tol_idx)-1, col_idx(1)+col_local(tol_idx)-1);
        peak_group(peak_neighbor) = group_num;
        peak_pool(peak_neighbor) = 0;
        
    end
    
    % Update list of current layer group members and output to preallocated
    % cell array
    layer_i = matrix_idx(peak_group==group_num);
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
layers = layers(cellfun(@(x) length(x) > 3, layers));

end