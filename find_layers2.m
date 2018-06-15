% Same as find_layer.m except that this function attempts to predict the
% next layer position based on previous values, rather than using the
% last known position

function [peak_group, layers] = find_layers2(peaks_raw, peak_width, core_res, horz_res)

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
    
    
    while search_R == true
        
        % Get nearest 10 group members
        group_idx = layer_i(end-min([19 length(layer_i)])+1:end);
        
        % Get row, col, and magnitude of current members of layer group
        [row_i, col_i] = ind2sub(size(peaks_raw), group_idx);
        mag_i = median(peaks_raw(group_idx));
        width_i = median(peak_width(group_idx));
        
%         if length(group_idx) <= 4
%             mag_var = 1;
%             
%         else
%             mag_var = var(peaks_raw(group_idx));
%         end
        
        switch length(group_idx) >= 10
            case false
                [row_n, ~] = ind2sub(size(peaks_raw), peak_n);
%                 col_n = max(col_i) + 1;
                col_n = max(col_i);
%                 row_var = var(row_i- mean(row_i));
            case true
                col_n = max(col_i);
%                 yy = smooth(col_i, row_i, 'rlowess');
%                 row_n = round(yy(end));
                
%                 Ri_mean = movmean(row_i, 5);
% %                 pp = csaps(col_i, Ri_mean, 0.25);
                pp = csaps(col_i, row_i);
%                 pp_val = fnval(pp, col_i);
                row_n = round(fnval(fnxtr(pp), col_n));
                
%                 row_var = var(row_i - pp_val);
        end
        
        
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
        
        dist_n = sqrt(((data_local(:,1)-row_n)/(0.50*width_i)).^2 + ...
            0.25*(data_local(:,2)-col_n).^2 + ((mag_local-mag_i)).^2);

        
        
%         if row_i <= 5
%             data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1)];
%             
%             dist_n = sqrt(((data_local(:,1)-row_n)/(0.50*width_i)).^2 + ...
%                 (data_local(:,2)-col_n).^2 + ((mag_local-mag_i)).^2);
%             
%         else
%             data_local = [(row_idx(1)+row_local-1)/(0.5*width_n) ...
%                 (col_idx(1)+col_local-1) mag_local];
%             data_i = [row_i col_i mag_i];
%             data_EX = repmat([row_n col_n mag_n], length(local_idx), 1);
%             sigma = cov(data_i);
%             
%             dist_n = sqrt(diag((data_local - data_EX)*inv(sigma)*...
%                 (data_local - data_EX)'));
%             %             dist_n = pdist2(data_local, data_EX(1,:), 'mahalanobis', sigma);
%             
%         end
        
        
        
        % Set distance threshold (based on p <= 0.01 for n-1 df on
        % Chi-squared distribution)
        %         threshold = 5.991;
        threshold = 6;
        [min_dist, dist_idx] = min(dist_n);
        
        if min_dist <= threshold
            
            peak_neighbor = sub2ind(size(peaks_raw), ...
                row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
            
            peak_group(peak_neighbor) = group_num;
            peak_pool(peak_neighbor) = 0;
            
            % Update list of current group members
            peak_n = peak_neighbor;
            layer_i = matrix_idx(peak_group==group_num);
            
            [~, col_check] = ind2sub(size(peaks_raw), peak_n);
            
            % If layer has reached the right edge of radargram, end search
            % to right of layer
            if col_check >= size(peaks_raw, 2)
                search_R = false;
            end
            
        else
            search_R = false;
        end
        
    end
    
    
    
    peak_n = peak_max;
    search_L = true;
    
    while search_L == true
        
        % Get nearest 10 group members
        group_idx = layer_i(1:min([19 length(layer_i)]));
        
        % Get row, col, and magnitude of current members of layer group
        [row_i, col_i] = ind2sub(size(peaks_raw), group_idx);
        mag_i = median(peaks_raw(group_idx));
        width_i = median(peak_width(group_idx));
        
        
%         if length(group_idx) <= 4
%             mag_var = 1;
%             
%         else
%             mag_var = var(peaks_raw(group_idx));
%         end
        
        switch length(group_idx) >= 10
            case false
                [row_n, ~] = ind2sub(size(peaks_raw), peak_n);
                col_n = min(col_i);
%                 col_n = min(col_i) - 1;
%                 row_var = var(row_i - mean(row_i));
                
            case true
                col_n = min(col_i);
%                 col_n = min(col_i) - 1;

%                 yy = smooth(col_i, row_i, 'rlowess');
%                 row_n = round(yy(1));
                
                pp = csaps(col_i, row_i);
%                 pp_val = fnval(pp, col_i);
                row_n = round(fnval(fnxtr(pp), col_n));
                
%                 row_var = var(row_i - pp_val);
%                 Ri_mean = movmean(row_i, 5);
% %                 pp = csaps(col_i, Ri_mean, 0.25);
        end

        
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
        
        % MAY NEED TO CHANGE TO LOCAL GROUP MEMBERS, NOT ALL GROUP MEMBERS
        %         data_i = [row_i col_i mag_i];
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1) mag_local];
        %         data_EX = repmat([row_pred col_pred mag_n], length(local_idx), 1);
        %         sigma = cov(data_i);
        
        dist_n = sqrt(((data_local(:,1)-row_n)/(0.5*width_i)).^2 + ...
            0.25*(data_local(:,2)-col_n).^2 + ((mag_local-mag_i)).^2);
        %         dist_n = sqrt(diag((data_local - data_EX)*inv(sigma)*...
        %             (data_local - data_EX)'));
        %         dist_n = pdist2(data_local, data_EX(1,:), 'mahalanobis', sigma);
        
        % Select the nearest neighbor to the estimated layer position
        %         [min_dist, dist_idx] = min(dist_n);
        
%         threshold = 5.991;
        threshold = 6;
        [min_dist, dist_idx] = min(dist_n);
        
        if min_dist <= threshold
            
            peak_neighbor = sub2ind(size(peaks_raw), ...
                row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
            
            peak_group(peak_neighbor) = group_num;
            peak_pool(peak_neighbor) = 0;
            
            % Update list of current group members
            peak_n = peak_neighbor;
            layer_i = matrix_idx(peak_group==group_num);
            
            [~, col_check] = ind2sub(size(peaks_raw), peak_n);
            
            % If layer has reached the left edge of radargram, end search
            % to left of layer
            if col_check <= 1
                search_L = false;
            end
            
        else
            search_L = false;
        end
    end
    

    
    % Find subscript indices of all members of ith layer
    [row, col] = ind2sub(size(peaks_raw), layer_i);
    
    % Smooth layer using a moving mean of row positions
    row_mean = round(movmean(row, 5));
    
    % For loop to search for additional peaks near current layer
    for j = 1:length(layer_i)
        
        % Find 10 nearest layer members on either side of jth member
        group_idx = sub2ind(size(peaks_raw), ...
            row(max([1 j-10]):min([length(row) j+10])), ...
            col(max([1 j-10]):min([length(col) j+10])));
        
%         row_group = row(max([1 j-10]):min([length(row) j+10]));
%         row_var = var(row_group - ...
%             row_mean(max([1 j-10]):min([length(row) j+10])));
        
        % Find local mean peak magnitude and width from the nearby layer
        % members
        mag_j = median(peaks_raw(group_idx));
        width_j = median(peak_width(group_idx));
        
        % Define local search window as 0.50 m above and below jth member
        % and 100 m to left and right of jth member
        row_idx = [max([row_mean(j)-round(0.50/core_res) 1]) ...
            min([row_mean(j)+round(0.50/core_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col(j)-round(150/horz_res)]) ...
            min([size(peaks_raw, 2) col(j)+round(150/horz_res)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1)];
        
        % Calculate distance between jth layer member and peaks within the
        % local search window
        dist_j = sqrt(((data_local(:,1)-row_mean(j))/(0.5*width_j)).^2 + ...
            (data_local(:,2)-col(j)).^2 + (1*(mag_local - mag_j)).^2);
%         dist_j = sqrt(((data_local(:,1)-row_mean(j))/(0.5*width_j+row_var)).^2 + ...
%             ((data_local(:,2)-col(j))/(0.5*(col_idx(2)-col_idx(1)))).^2);
        
        % Set distance threshold based on peak width and error bin size
        threshold = 3;
        
        % Assign all peaks within tolerance to the current layer group and
        % remove those peaks from the pool of available peaks to search
        tol_idx = dist_j <= threshold;
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(tol_idx)-1, col_idx(1)+col_local(tol_idx)-1);
        peak_group(peak_neighbor) = group_num;
        peak_pool(peak_neighbor) = 0;
        
    end
    
    
    layer_i = matrix_idx(peak_group==group_num);
    
    
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

% Remove empty cells and layers shorter than 250 m
layers = layers(cellfun(@(x) length(x) > 10, layers));

end