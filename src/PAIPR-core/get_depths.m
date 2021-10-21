% Function to find the depth of each assigned year layer

function [yr_DepthTables] = get_depths(radar, bin_size)

% Determine number of elements in each bin based on distance spacing and
% requested bin_size (bin_size in meters)
dist_interval = mean(diff(radar.dist));
stack_step = round(bin_size / dist_interval);

% Determine starting and stopping indices for different aggregated bins
start_idx = 1:stack_step:length(radar.SMB);
end_idx = [start_idx(2:end)-1 length(radar.SMB)];

% Preallocate cell array for annual layer depth tables
yr_DepthTables = cell(1, length(start_idx));


for i=1:length(yr_DepthTables)
    
    % Extract ages with current bin and reshape to 2D
    ages_i = radar.ages(:,start_idx(i):end_idx(i),:);
    ages_i = reshape(ages_i, size(ages_i,1), []);

    % Get min/max ages in slice
    yr_first = floor(max(max(ages_i)));
    yr_last = ceil(min(min(ages_i)));
    yr_list = (yr_first:-1:yr_last)';
    
    % Logical of integer year indices
    yr_logic = logical(diff(floor(ages_i)));
    % Preallocate array for depth of annual IRHs
    yr_depths = NaN(length(yr_list), size(ages_i,2));
    for j=1:size(yr_logic,2)
        
        % Get depths of annual IRHs for each simulation
        yr_idx = find(yr_logic(:,j));
        yr_depths(1:length(yr_idx),j) = radar.depth(yr_idx);
    end
    
    % Remove years with more than 25% NaN values
    nan_idx = sum(isnan(yr_depths),2)/size(yr_depths,2);
    yr_list = yr_list(nan_idx <= 0.25);
    yr_depths = yr_depths(nan_idx<=0.25,:);
    
    % Determine average bin values for remaining variables of interest
    yr_len = length(yr_list);
    time = mean(radar.collect_time(start_idx(i):end_idx(i)));
    lat = mean(radar.Lat(start_idx(i):end_idx(i)));
    lon = mean(radar.Lon(start_idx(i):end_idx(i)));
    elev = mean(radar.elev(start_idx(i):end_idx(i)));
    
    % Output years, median depth, and std. err for each bin
    yr_DepthTables{i} = table(...
        repmat(time, yr_len,1), repmat(lat, yr_len,1), ...
        repmat(lon, yr_len,1), repmat(elev, yr_len,1), ...
        yr_list, median(yr_depths,2, 'omitnan'), ...
        std(yr_depths,1,2,'omitnan'), ...
        sum(~isnan(yr_depths),2),...
        'VariableNames', {'collect_time', 'Lat', 'Lon', 'elev', ...
        'IHR_year', 'IHR_depth', 'IHR_std', 'IHR_cnt'});
    
    
end





end