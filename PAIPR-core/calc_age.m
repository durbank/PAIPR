function [radar] = calc_age(radar_struct, cores, Ndraw)

% Convert to depth
[radar] = radar_depth(radar_struct, cores);

% Find the mean response with depth in the radar data attributes across a
% given horizontal resolution (in meters)
horz_res = 25;
[radar] = radar_stack(radar, horz_res);

%% Signal-noise processing

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

% Define the vertical resolution of the core data
core_res = 0.02;

% Define the cutoff depth for radar traces and find index of crossover
% depth
% cutoff = 25;
cutoff = 30;
depth_bott = floor(min([min(radar.depth(end,:)) cutoff]));

% Trim radar traces to cutoff depth and interpolate data to vertical scale
% of the firn cores
radarZ_interp = zeros(depth_bott/core_res+1, size(radar.data_stack, 2));
for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:core_res:radar.depth(end,i));
    radarZ_i = interp1(radar.depth(:,i), radar_Z(:,i), depth_interp, 'pchip');
    radarZ_interp(:,i) = radarZ_i(1:size(radarZ_interp, 1));
end

% % If manual layer picks are present, transform them to same depth and
% % vertical scale as the interpolated radar data
% if isfield(radar, 'man_layers')
%     man_interp = zeros(size(radarZ_interp));
%     for i = 1:size(radar.data_stack, 2)
% %         depth_interp = (0:core_res:radar.depth(end,i))';
%         layer_idx = logical(radar.man_layers(:,i));
%         layer_num = radar.man_layers(layer_idx,i);
%         man_depth_i = radar.depth(layer_idx,i);
%         depth_idx = round(man_depth_i/core_res) + 1;
%         man_interp(depth_idx(man_depth_i<=cutoff),i) = layer_num(man_depth_i<=cutoff);
%     end
%     radar.man_layers = man_interp;
% end

% Assign structure output depth to interpolated depths
radar.depth = (0:core_res:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~18 cm)
radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);

% Clear unnecessary variables
clearvars -except file cores Ndraw radar horz_res core_res

%%

% Iterative radon transforms
[IM_gradients] = radar_gradient(radar, core_res, horz_res);

% Find radar peaks in echogram
[peaks_raw, peak_width] = radar_peaks(radar, core_res);

% Find continuous layers within radargram based on peaks and layer stream
% field
[peaks, group_num, layers] = radar_trace(peaks_raw, peak_width, ...
    IM_gradients, core_res, horz_res);


% %% Clean layer picks
% 
% % Preallocate arrays for the matrix indices of members of each layer
% layers_idx = cell(1,length(layers));
% peaks = zeros(size(peaks_raw));
% 
% % For loop to coerce layers to have one row position for each trace
% for i = 1:length(layers_idx)
%     
%     % Find matrix indices of all members of ith layer
%     layer_i = layers{i};
%     
%     % Find row and col indices of members of ith layer
%     [row, col] = ind2sub(size(radar.data_smooth), layer_i);
%     mag = peaks_raw(layer_i);
%     
%     % Interpolate data to all column positions within the range of the
%     % layer
%     col_interp = min(col):max(col);
%     
%     % Interpolate row positions using a cubic smoothing spline
%     row_interp = round(fnval(csaps(col, row), col_interp));
%     row_interp(row_interp < 1) = 1;
%     row_interp(row_interp > size(peaks,1)) = size(peaks,1);
%     
%     % Interpolate peak prominence magnitudes to all columns in range using
%     % a cubic smoothing spline
%     mag_interp = csaps(col, mag, 1/length(col_interp), col_interp);
%     
%     % Assign interpolated layer to output
%     layer_interp = sub2ind(size(peaks), row_interp, col_interp);
%     peaks(layer_interp) = mag_interp;
%     layers_idx{i} = layer_interp';
% end
% 
% % Create matrix of layer group assignments
% group_num = zeros(size(peaks));
% for i = 1:length(layers_idx)
%     group_num(layers_idx{i}) = i;
% end

% Output layer arrays to radar structure
radar.peaks = peaks;
radar.layers = layers;
radar.groups = group_num;

% Calculate age-depth profile distributions for each echogram trace
r = -3.06e-4;
k = 2.94;
[radar] = radar_age(radar, r, k, Ndraw);

% Clip depth-related variables to final cutoff depth
cutoff = 25;
cut_idx = min([(round(cutoff/core_res)+1) length(radar.depth)]);
radar_new = struct('collect_date', radar.collect_date, 'Easting', ...
    radar.Easting, 'Northing', radar.Northing, 'dist', radar.dist, ...
    'depth', radar.depth(1:cut_idx), 'rho_coeff', radar.rho_coeff, ...
    'rho_var', radar.rho_var, 'data_smooth', radar.data_smooth(1:cut_idx,:),...
    'peaks', radar.peaks(1:cut_idx,:), 'groups', radar.groups(1:cut_idx,:),...
    'likelihood', radar.likelihood(1:cut_idx,:), ...
    'ages', radar.ages(1:cut_idx,:,:));
if isfield(radar, 'elev')
    radar_new.elev = radar.elev;
end
radar = radar_new;

% Find layer member indices based on new clipped record
layers = cell(1, max(radar.groups(:)));
for j = 1:length(layers)
    layers{j} = find(radar.groups == j);
end
radar.layers = layers(~cellfun(@isempty, layers));

% Eventually will move the gamma-fitting process to here? In this way, we
% can significantly reduce the necessary storage space once we increase MC
% simulations

end