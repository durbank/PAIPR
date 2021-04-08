% Function to combine SMB MC results from nearby traces (within a given
% distance) into a single location distribution

function [accum_table] = accum_Gauss(radar, bin_size)

% radar = load('~/scratch/echo/radar_out1.mat');

dist_interval = mean(diff(radar.dist));
stack_step = round(bin_size / dist_interval);

start_idx = 1:stack_step:length(radar.SMB);
end_idx = [start_idx(2:end)-1 length(radar.SMB)];



% Perform gamma fitting for each iteration (saves memory load)
% For development, will do mean and std instead
accum_table = cell(1,length(start_idx));
% SMB_bar = cell(1,length(start_idx));
% SMB_std = SMB_bar;
% SMB_yr = SMB_bar;
% time = NaT(1,length(start_idx));
% lat = zeros(1,length(start_idx));
% lon = lat;
% elev = lat;


for i=1:length(start_idx)
    
    SMB_i = radar.SMB(start_idx(i):end_idx(i));
    yr_i = radar.SMB_yr(start_idx(i):end_idx(i));
    
    
    cell_len = cellfun(@length, yr_i);
    [~,max_idx] = max(cell_len);
    SMB_mat = nan(max(cell_len), length(SMB_i),size(SMB_i{1},2));
    for j=1:length(SMB_i)
        SMB_mat(1:cell_len(j),j,:) = SMB_i{j};
    end
    
    SMB_mat = reshape(...
        SMB_mat, size(SMB_mat,1), size(SMB_mat,2)*size(SMB_mat,3));
    nan_sum = sum(isnan(SMB_mat),2) / size(SMB_mat,2);
    SMB_mat = SMB_mat(nan_sum <= 0.33,:);
    
%     SMB_bar{i} = nanmean(SMB_mat,2);
%     SMB_std{i} = nanstd(SMB_mat,[],2);
%     SMB_yr{i} = yr_i{max_idx}(nan_sum <= 0.33);
    
    yr_len = size(SMB_mat,1);
    time = mean(radar.collect_time(start_idx(i):end_idx(i)));
    lat = mean(radar.Lat(start_idx(i):end_idx(i)));
    lon = mean(radar.Lon(start_idx(i):end_idx(i)));
    elev = mean(radar.elev(start_idx(i):end_idx(i)));
    
    
    accum_table{i} = table(repmat(time, yr_len,1), repmat(lat, yr_len,1), ...
        repmat(lon, yr_len,1), repmat(elev, yr_len,1), ...
        yr_i{max_idx}(nan_sum <= 0.33), nanmean(SMB_mat,2), ...
        nanstd(SMB_mat,[],2), 'VariableNames', ...
        {'collect_time', 'Lat', 'Lon', 'elev', ...
        'Year', 'accum_mu', 'accum_std'}); 
end
