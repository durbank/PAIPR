function [radar] = calc_age2(radar, r, k, Ndraw)
%STREAM_SUM Summary of this function goes here
%   Detailed explanation goes here

%% Signal-noise processing

% Stationarize the radar response by differencing traces with a smoothing
% spline
s = zeros(size(radar.data_stack));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth(:,i), radar.data_stack(:,i), ...
        0.95, radar.depth(:,i));
end
radar_stat = radar.data_stack - s;

% Define the vertical resolution of the core data and horizontal resolution
% of the radar data
vert_res = 0.02;
horz_res = round(mean(diff(radar.dist)));

% Define the cutoff depth for radar traces and find index of crossover
% depth
cutoff = 30;
depth_bott = floor(min([min(radar.depth(end,:)) cutoff]));

% Trim radar traces to cutoff depth and interpolate data to vertical scale
% of the firn cores
radar_interp = zeros(depth_bott/vert_res+1, size(radar.data_stack, 2));
for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:vert_res:radar.depth(end,i));
    radar_i = interp1(radar.depth(:,i), radar_stat(:,i), ...
        depth_interp, 'pchip');
    radar_interp(:,i) = radar_i(1:size(radar_interp, 1));
end

% Assign structure output depth to interpolated depths
radar.depth = (0:vert_res:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~18 cm)
radar.data_smooth = sgolayfilt(radar_interp, 3, 9);




% Iterative radon transforms
[IM_gradients] = radar_gradient(radar, vert_res, horz_res);

% 
[stream_val, XY_stream] = stream_sum(radar, IM_gradients);


% 
[~, peak_idx, ~, peak_prom] = findpeaks(stream_val, ...
        'MinPeakProminence', std(stream_val)/10);

%%

% Find the mean peak prominence (used to scale the prominence-distance
% results)
peak_w = 1/std(radar.data_smooth(:));
dist_w = 1/(size(radar.data_smooth,2)*horz_res);
% peak_w = 1;
% dist_w = 1;


layer_DB = zeros(size(radar.data_smooth));
for i=1:length(peak_idx)
    
    layer_DB(XY_stream{peak_idx(i)}) = peak_w*dist_w*peak_prom(i);
end

% Preallocate arrays for layer likelihoods and anges
ages = zeros([size(radar.data_smooth) Ndraw]);
radar.likelihood = zeros(size(radar.data_smooth));
err_out = [];
for i=1:size(layer_DB,2)
    
    % Define surface age and the year associated with the first pick of the
    % algorithm
    yr_vec = datevec(radar.collect_time(i));
    yr_pick1 = yr_vec(1);
    age_top = yr_vec(1) + (30*yr_vec(2)+yr_vec(3))/365;
    
    % Get layer DB values and depths for layers in ith trace
    layer_i = layer_DB(:,i);
    layer_idx = layer_i>0;
    layer_i = layer_i(layer_idx);
    depths_i = radar.depth(layer_idx);
    
    % 
    likelihood = 1./(1+exp(r*layer_i + k));
    radar.likelihood(layer_idx,i) = likelihood;
    
    % Assign MC simulation annual layer presence based on layer likelihood
    % values
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end
    
    for j = 1:Ndraw
        depths_j = [0; depths_i(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        try
            ages(:,i,j) = interp1(depths_j, yrs_j, radar.depth, ...
                'linear', 'extrap');
        catch
            sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
            err_out = [err_out j];
        end
    end
    if ~isempty(err_out)
        ages(:,i,err_out) = repmat(sum(squeeze(ages(:,i,:)), 2)./...
            sum(squeeze(ages(:,i,:))~=0, 2), 1, length(err_out));
    end
    err_out = [];
end

radar.ages = ages;

%%

% Clip depth-related variables to final cutoff depth
cutoff = 25;
cut_idx = min([(round(cutoff/vert_res)+1) length(radar.depth)]);
radar_new = struct('collect_time', radar.collect_time, 'Easting', ...
    radar.Easting, 'Northing', radar.Northing, 'dist', radar.dist, ...
    'depth', radar.depth(1:cut_idx), 'data_smooth', ...
    radar.data_smooth(1:cut_idx,:), 'DB', layer_DB(1:cut_idx,:), ...
    'likelihood', radar.likelihood(1:cut_idx,:), ...
    'ages', radar.ages(1:cut_idx,:,:));

if isfield(radar, 'elev')
    radar_new.elev = radar.elev;
end
radar = radar_new;


end