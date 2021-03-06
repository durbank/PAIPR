% Second refined method for picking continuous layers. Methods used are
% similar to the nearest neighbor search used in find_layers.m, but this
% method uses streamlines generated from a layer's depth-slope to estimate
% the position of the nearest neighbors rather than relying on the position
% of the last member found

function [peaks, group_num, layers] = radar_trace(peaks_raw, ...
    peak_width, grad_matrix, vert_res, horz_res)

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

% Calculate 
xstart_idx = round(size(peaks_raw,2)/2);
[stream_val, ~] = stream_sum(peaks_raw, grad_matrix, ...
    xstart_idx, horz_res);
min_Prom = quantile(stream_val, 0.05);

% While loop that iterates on each accumulation layer
while search_new == true
    
    
    % 
    [stream_val, XY_streams] = stream_sum(peak_pool, grad_matrix, ...
    xstart_idx, horz_res);
    [~, peak_idx, ~, peak_prom] = findpeaks(stream_val, ...
        'MinPeakProminence', min_Prom, ...
        'MinPeakDistance', round(0.08/vert_res));
    
    if isempty(peak_prom)
        break
    end
    
    [max_val, max_idx] = max(peak_prom);
    stream_idx = XY_streams{peak_idx(max_idx)};
    stream_peaks = peak_pool(stream_idx);
    [~, peak_max] = max(stream_peaks);
    
    if max_val <= min_Prom
        search_new = false;
    end
    
    % Assign max peak to current group and remove from future peak searches
    peak_group(stream_idx(peak_max)) = group_num;
    peak_pool(stream_idx(peak_max)) = 0;
    
    
%     % Convolve remaining peaks with a 10x10 matrix and find brightest
%     % convolved pixel (2D smoothing of matrix in order to find brightest
%     % local region)
%     A = peak_pool;
%     B = ones(10)/10^2;
%     C = conv2(A,B,'same');
%     [~,C_max] = max(C(:));
%     [r_C, c_C] = ind2sub(size(C), C_max);
%     
%     % Find local area around max convolved peak
%     r_idx = [max([1 r_C-5]) min([size(peaks_raw,1) r_C+5])];
%     c_idx = [max([1 c_C-5]) min([size(peaks_raw,2) c_C+5])];
%     
%     % Find index of actual brightest peak within local search area
%     peak_local = peak_pool(r_idx(1):r_idx(2),c_idx(1):c_idx(2));
%     [~, local_max] = max(peak_local(:));
%     [r_max, c_max] = ind2sub(size(peak_local), local_max);
%     peak_max = sub2ind(size(peak_pool), ...
%         r_idx(1)+r_max-1, c_idx(1)+c_max-1);
%     
%     % Assign max peak to current group and remove from future peak searches
%     peak_group(peak_max) = group_num;
%     peak_pool(peak_max) = 0;
    
    %% Righthand search for layer members
    
    % Initialize values for group search within layer to right of max peak
    layer_i = matrix_idx(peak_group==group_num);
    search_R = true;
    
    threshold = 3;
    
    while search_R == true
        
        % Get nearest 20 group members
        group_idx = layer_i(end-min([20 length(layer_i)])+1:end);
        
        % Get median magnitude and median peak width of 10 nearest group
        % members
        mag_i = mean(peaks_raw(group_idx));
        Hwidth_i = 0.5*mean(peak_width(group_idx));
        
        % Get subscripts of farthest right 20 group members
        [row_i, col_i] = ind2sub(size(peaks_raw), group_idx(end - ...
            min([20 length(layer_i)])+1:end));
        
        % Get col and row subscripts of farthest right position of layer_i
        col_n = max(col_i);
%         row_n = round(mean(row_i));
        row_n = row_i(end);
        
        % Define local  search window as 250 m to right of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_n-round(0.50/vert_res) 1]) ...
            min([row_n+round(0.50/vert_res) size(peaks_raw, 1)])];
        col_idx = [col_n min([col_n+round(250/horz_res) ...
            size(peaks_raw, 2)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(...
            row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1)];
        
        % Estimate layer streams within the local search window for nearby
        % members of layer_i
        stream_n = stream2(ones(size(peaks_raw)), grad_matrix, ...
            col_i, row_i, [0.10, 1000]);
        cols_stream = col_idx(1):col_idx(2);
        rows_stream = zeros(length(stream_n), size(peak_local,2));
        for k = 1:length(stream_n)
            for kk = 1:length(cols_stream)
               kk_start = find(stream_n{k}(:,1)>=cols_stream(kk), 1);
               kk_end = find(stream_n{k}(:,1)<=...
                   cols_stream(kk)+1, 1, 'last');
               rows_stream(k,kk) = mean(stream_n{k}(kk_start:kk_end,2));
            end
        end
        
        % Estimate the average layer stream positions for layer_i within
        % the local search window
        data_stream = [mean(rows_stream, 1)' cols_stream'];
        
        % Calculate the weights to use on the indiviudual dissimilarity
        % components
        w_x = (1/std(col_local))^2;
        w_y = (1/(Hwidth_i))^2;
        if length(group_idx) <= 5
            mag_i_std = 0.10*mean(peaks_raw(group_idx));
        else
            mag_i_std = std(peaks_raw(group_idx));
        end
        w_m = (1/mag_i_std)^2;

        % Calculate the standardized dissmilarity values for each local
        % peak
        dist_n = sqrt(w_x*(data_local(:,2)-(data_stream(:,2))').^2 + ...
            w_y*(data_local(:,1)-(data_stream(:,1))').^2 + ...
            w_m*repmat(mag_local-mag_i, 1, size(data_stream,1)).^2);
        
        % Find all peaks within tolerance to the layer streamline
        dist_idx = min(dist_n, [], 2) <= threshold;
        
        % Get matrix index of nearest neighbors
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
        
        % Assign nearest neighbors to layer group, and remove peaks from
        % future searches
        peak_group(peak_neighbor) = group_num;
        peak_pool(peak_neighbor) = 0;
        
        % Update list of current group members
        layer_i = matrix_idx(peak_group==group_num);
        
        % Check if current layer has reached the right edge of the radar
        % gram
        [~, col_check] = ind2sub(size(peaks_raw), layer_i(end));
        
        % If layer has no nearest neighbors within tolerance or has reached
        % the right edge of radargram, end search to right of layer
        if isempty(peak_neighbor) || col_check >= size(peaks_raw, 2)
            search_R = false;
        end
        
    end

    %% Lefthand search for layer members
    
    % Intialize peak search to the left
    search_L = true;
    
    while search_L == true
        
        % Get nearest 20 group members
        group_idx = layer_i(1: min([20 length(layer_i)]));
        
        % Get median magnitude and median peak width of 10 nearest group
        % members
        mag_i = mean(peaks_raw(group_idx));
        Hwidth_i = 0.5*mean(peak_width(group_idx));
        
        % Get subscripts of farthest left 20 group members
        [row_i, col_i] = ind2sub(size(peaks_raw), ...
            group_idx(1:min([20 length(layer_i)])));
        
        % Get col and row subscripts of farthest left position of layer_i
        col_n = min(col_i);
%         row_n = round(mean(row_i));
        row_n = row_i(1);
        
        % Define local  search window as 250 m to right of last known group
        % member, and 0.50 m above and below estimated row position for the
        % projected nearest neighbor
        row_idx = [max([row_n-round(0.50/vert_res) 1]) ...
            min([row_n+round(0.50/vert_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col_n-round(250/horz_res)]) col_n];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1)];
        
        % Estimate layer streams (leftward moving) within the local search 
        % window for nearby members of layer_i
        stream_n = stream2(-ones(size(peaks_raw)), -grad_matrix, ...
            col_i, row_i, [0.10, 1000]);
        cols_stream = col_idx(1):col_idx(2);
        rows_stream = zeros(length(stream_n), size(peak_local,2));
        for k = 1:length(stream_n)
            for kk = 1:length(cols_stream)
               kk_start = find(stream_n{k}(:,1)<=cols_stream(kk)+1, 1);
               kk_end = find(stream_n{k}(:,1)>=cols_stream(kk), 1, 'last');
               rows_stream(k,kk) = mean(stream_n{k}(kk_start:kk_end,2));
            end
        end
        
        % Estimate the average layer stream positions for layer_i within
        % the local search window
        data_stream = [mean(rows_stream, 1)' cols_stream'];
        
        % Calculate the weights to use on the indiviudual dissimilarity
        % components
        w_x = (1/std(col_local))^2;
        w_y = (1/(Hwidth_i))^2;
        if length(group_idx) <= 5
            mag_i_std = 0.10*mean(peaks_raw(group_idx));
        else
            mag_i_std = std(peaks_raw(group_idx));
        end
        w_m = (1/mag_i_std)^2;

        % Calculate the standardized dissmilarity values for each local
        % peak
        dist_n = sqrt(w_x*(data_local(:,2)-(data_stream(:,2))').^2 + ...
            w_y*(data_local(:,1)-(data_stream(:,1))').^2 + ...
            w_m*repmat(mag_local-mag_i, 1, size(data_stream,1)).^2);
        
        % Find all peaks within tolerance to the layer streamline
        dist_idx = min(dist_n, [], 2) <= threshold;
        
        % Get matrix index of nearest neighbors
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(dist_idx)-1, col_idx(1)+col_local(dist_idx)-1);
        
        % Assign nearest neighbors to layer group, and remove peaks from
        % future searches
        peak_group(peak_neighbor) = group_num;
        peak_pool(peak_neighbor) = 0;
        
        % Update list of current group members
        layer_i = matrix_idx(peak_group==group_num);
        
        % Check if current layer has reached the left edge of the radargram
        [~, col_check] = ind2sub(size(peaks_raw), layer_i(1));
        
        % If layer has no nearest neighbors within tolerance or has reached
        % the left edge of radargram, end search to left of layer
        if isempty(peak_neighbor) || col_check <= 1
            search_L = false;
        end
        
    end
    
    %% Add additional nearby peaks to layer
    % (Those that were not nearest neighbors but are still sufficiently
    % close to be included in the current layer)

    % Find subscript indices of all members of ith layer
    [row, col] = ind2sub(size(peaks_raw), layer_i);
    
    % Smooth layer using a moving mean of row positions
    row_mean = round(movmean(row, 5));
    
    % For loop to search for additional peaks near current layer
    for j = 1:length(layer_i)
        
        % Find 10 nearest layer members on either side of jth member
        group_idx = sub2ind(size(peaks_raw), ...
            row(max([1 j-5]):min([length(row) j+5])), ...
            col(max([1 j-5]):min([length(col) j+5])));
        
        % Find local mean peak magnitude and width from the nearby layer
        % members
        mag_j = mean(peaks_raw(group_idx));
        Hwidth_j = 0.5*mean(peak_width(group_idx));
        
        % Define local search window as 0.50 m above and below jth member
        % and 250 m to left and right of jth member
        row_idx = [max([row_mean(j)-round(0.50/vert_res) 1]) ...
            min([row_mean(j)+round(0.50/vert_res) size(peaks_raw, 1)])];
        col_idx = [max([1 col(j)-round(250/horz_res)]) ...
            min([size(peaks_raw, 2) col(j)+round(150/horz_res)])];
        
        % Get matrix of peaks from avaiable pool within the local search
        % area, along with their indices, magnitudes, and subscript indices
        peak_local = peak_pool(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        local_idx = find(peak_local);
        mag_local = peak_local(local_idx);
        [row_local, col_local] = ind2sub(size(peak_local), local_idx);
        data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1)];
        
        % Calculate the weights to use on the indiviudual dissimilarity
        % components
        w_x = (1/std(col_local))^2;
        w_y = (1/(Hwidth_j))^2;
        if length(group_idx) <= 5
            mag_j_std = 0.10*mean(peaks_raw(group_idx));
        else
            mag_j_std = std(peaks_raw(group_idx));
        end
        w_m = (1/mag_j_std)^2;
        w_m = 0;    % Set to zero so all peaks within distance are added
                    % regardless of magnitude differences

        % Calculate the standardized dissmilarity values for each local
        % peak
        dist_j = sqrt(w_x*(data_local(:,2)-col(j)).^2 + ...
            w_y*(data_local(:,1)-row_mean(j)).^2 + ...
            w_m*(mag_local-mag_j).^2);
        
        % Set distance threshold for secondary recovery search
        threshold = 1.0;
        
        % Assign all peaks within tolerance to the current layer group and
        % remove those peaks from the pool of available peaks to search
        tol_idx = dist_j <= threshold;
        peak_neighbor = sub2ind(size(peaks_raw), ...
            row_idx(1)+row_local(tol_idx)-1, col_idx(1)+col_local(tol_idx)-1);
        peak_group(peak_neighbor) = group_num;
        peak_pool(peak_neighbor) = 0;
        
    end
    
    % Update list of current group members
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


%% Clean final layer picks

% Remove empty cells and layers shorter than 500 m
layers = layers(cellfun(@(x) length(x) >= round(500/horz_res), layers));

% Preallocate arrays for the matrix indices of members of each layer
layers_idx = cell(1,length(layers));
peaks = zeros(size(peaks_raw));

% For loop to coerce layers to have one row position for each trace
for i = 1:length(layers_idx)
    
    % Find matrix indices of all members of ith layer
    layer_i = layers{i};
    
    % Find row and col indices of members of ith layer
    [row, col] = ind2sub(size(peaks), layer_i);
    mag = peaks_raw(layer_i);
    
    % Interpolate data to all column positions within the range of the
    % layer
    col_interp = min(col):max(col);
    
    % Interpolate row positions using a cubic smoothing spline
    row_interp = round(fnval(csaps(col, row), col_interp));
    row_interp(row_interp < 1) = 1;
    row_interp(row_interp > size(peaks,1)) = size(peaks,1);
    
    % Interpolate peak prominence magnitudes to all columns in range using
    % a cubic smoothing spline
    mag_interp = csaps(col, mag, 1/length(col_interp), col_interp);
    
    % Assign interpolated layer to output
    layer_interp = sub2ind(size(peaks), row_interp, col_interp);
    peaks(layer_interp) = mag_interp;
    layers_idx{i} = layer_interp';
end

% Create matrix of layer group assignments
group_num = zeros(size(peaks));
for i = 1:length(layers_idx)
    group_num(layers_idx{i}) = i;
end

% Rename layer cell array for output
layers = layers_idx;

end