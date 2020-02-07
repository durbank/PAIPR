function [radar] = calc_age(radar, r, k, Ndraw)



%% Signal-noise processing

% Stationarize the radar response by differencing traces with a smoothing
% spline
s = zeros(size(radar.data_stack));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth(:,i), radar.data_stack(:,i), ...
        0.95, radar.depth(:,i));
end
radar_stat = radar.data_stack - s;

% Remove linear trend in variance (attentuation with depth) and convert to
% standardized values (z-score statistics)
radar_Z = zeros(size(radar_stat));
for i = 1:size(radar_stat, 2)
    data_i = radar_stat(:,i);
    % Frame length to define local variance
    half_frame = round(0.5*length(data_i)/5);
    var0 = movvar(data_i, 2*half_frame, 'EndPoints', 'discard');
    x = (half_frame:length(data_i)-half_frame)';
    % Linear trend in variance
    EQ = polyfit(x, var0, 1);
    x_mod = (1:length(data_i))';
    mod = polyval(EQ, x_mod);
    % Standardize variance-corrected data
    radar_Z(:,i) = data_i./sqrt(abs(mod));
end

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
radarZ_interp = zeros(depth_bott/vert_res+1, size(radar.data_stack, 2));
for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:vert_res:radar.depth(end,i));
    radarZ_i = interp1(radar.depth(:,i), radar_Z(:,i), ...
        depth_interp, 'pchip');
    radarZ_interp(:,i) = radarZ_i(1:size(radarZ_interp, 1));
end

% Assign structure output depth to interpolated depths
radar.depth = (0:vert_res:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~18 cm)
radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);

% Clear unnecessary variables
clearvars -except radar r k Ndraw horz_res vert_res
%%

% Iterative radon transforms
[IM_gradients] = radar_gradient(radar, vert_res, horz_res);

% Find radar peaks in echogram
[peaks_raw, peak_width] = radar_peaks(radar, vert_res);



[peaks, group_num, layers] = radar_trace(peaks_raw, peak_width, ...
    IM_gradients, vert_res, horz_res);

% Find max depth with full layer coverage
max_depth = zeros(1,size(group_num,2));
for i=1:length(max_depth)
    max_depth(i) = find(group_num(:,i), 1, 'last');
end
cut_idx = min(max_depth);

[rows,cols] = cellfun(@(x) ind2sub(size(peaks),x), layers, ...
    'UniformOutput', false);
[layers_new] = cellfun(@(x,y) sub2ind([cut_idx length(radar.dist)], ...
    x(x<=cut_idx), y(x<=cut_idx)), rows, cols, 'UniformOutput', false);
layers_new = layers_new(~cellfun(@isempty, layers_new));


radar.depth = radar.depth(1:cut_idx);
radar.data_smooth = radar.data_smooth(1:cut_idx,:);
radar.IM_grad = IM_gradients(1:cut_idx,:);
radar.peaks = peaks(1:cut_idx,:);
radar.layers = layers_new;
radar.groups = group_num(1:cut_idx,:);
radar = rmfield(radar, {'data_stack', 'time_trace'});

% Calculate age-depth profile distributions for each echogram trace
[radar] = radar_age(radar, r, k, Ndraw);

% % Clip depth-related variables to final cutoff depth
% cutoff = 25;
% cut_idx = min([(round(cutoff/vert_res)+1) length(radar.depth)]);
% radar_new = struct('collect_time', radar.collect_time, 'Easting', ...
%     radar.Easting, 'Northing', radar.Northing, 'dist', radar.dist, ...
%     'depth', radar.depth(1:cut_idx), 'data_smooth', ...
%     radar.data_smooth(1:cut_idx,:), 'peaks', radar.peaks(1:cut_idx,:), ...
%     'groups', radar.groups(1:cut_idx,:), 'likelihood', ...
%     radar.likelihood(1:cut_idx,:), 'ages', radar.ages(1:cut_idx,:,:));
% 
% if isfield(radar, 'elev')
%     radar_new.elev = radar.elev;
% end
% radar = radar_new;
% 
% % Find layer member indices based on new clipped record
% layers_new = cell(1, max(radar.groups(:)));
% for j = 1:length(layers_new)
%     layers_new{j} = find(radar.groups == j);
% end
% radar.layers = layers_new(~cellfun(@isempty, layers_new));

end