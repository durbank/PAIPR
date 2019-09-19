function [radar] = radar_RT(radar_struct, cores, Ndraw)

% Convert to depth
[radar] = radar_depth(radar_struct, cores);

% Find the mean response with depth in the radar data attributes across a
% given horizontal resolution (in meters)
horz_res = 25;
[radar] = radar_stack(radar, horz_res);

%%

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
[IM_gradients] = radar_radon(radar, core_res, horz_res);

% Find radar peaks in echogram
[peaks_raw, peak_width] = radar_peaks(radar);

% Find continuous layers within radargram based on peaks and layer stream
% field
[~, layers] = find_layers2(peaks_raw, peak_width, ...
    IM_gradients, core_res, horz_res);



%%

% Preallocate arrays for the matrix indices of members of each layer
layers_idx = cell(1,length(layers));
peaks = zeros(size(peaks_raw));

% For loop to coerce layers to have one row position for each trace
for i = 1:length(layers_idx)
    
    % Find matrix indices of all members of ith layer
    layer_i = layers{i};
    
    % Find row and col indices of members of ith layer
    [row, col] = ind2sub(size(radar.data_smooth), layer_i);
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

%%

% Calculate continuous layer distances for each layer (accounting for 
% lateral size of stacked radar trace bins)
layers_dist = cellfun(@(x) numel(x)*horz_res, layers_idx);

% Map layer prominence-distance values to the location within the radar
% matrix of the ith layer
layer_peaks = zeros(size(peaks));
for i = 1:length(layers_idx)
    layer_peaks(layers_idx{i}) = peaks(layers_idx{i}).*layers_dist(i);
end


% Output layer arrays to radar structure
radar.peaks = peaks;
radar.layers = layers_idx;
radar.groups = group_num;

%% Assign layer likelihood scores and estimate age-depth scales

% Define surface age and the year associated with the first pick of the 
% algorithm
age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);

% Preallocate arrays for layer likelihoods and anges
ages = zeros([size(radar.data_smooth) Ndraw]);
radar.likelihood = zeros(size(radar.data_smooth));
err_out = [];
for i = 1:size(layer_peaks, 2)
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = layer_peaks(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar.depth(peaks_idx);
    
    % Likelihood of layer representing a year based on a logistic function
    % with rate (r) calculated above
%     r = -2.4333e-4; % [-3.18e-4 -1.55e-4]
%     k = 4.4323;     % [3.25 4.8]
    r = -3.06e-4;
    k = 2.94;
    
    likelihood = 1./(1+exp(r*peaks_i + k));
    radar.likelihood(peaks_idx,i) = likelihood;
    
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
            ages(:,i,j) = interp1(depths_j, yrs_j, radar.depth, 'linear', 'extrap');
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
end