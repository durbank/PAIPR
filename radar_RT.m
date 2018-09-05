function [radar] = radar_RT(radar_file, cores, Ndraw)

radar_type = isstruct(radar_file);

switch radar_type
    case false
        % Conversion to depth
        [radar] = radar_depth(radar_file, cores);
        
    case true
        % Conversion to depth
        [radar] = OIB_depth(radar_file, cores);
end

% % Determine if there are data break points (large gaps in data that would
% % necessitate data processing over segments of the whole data)
% distance = [0 diff(radar.dist)];
% dist_idx = distance >= 500;
% data_col = 1:size(radar.data_out, 2);
% data_endpts = [1 data_col(dist_idx)-1 length(data_col)];

% Find the mean response with depth in the radar data attributes across a
% given horizontal resolution (in meters)
horz_res = 25;
[radar] = radar_stack(radar, horz_res);

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
    data = radar_stat(:,i);
    % Frame length to define local variance
    half_frame = round(0.5*length(data)/5); 
    var0 = movvar(data, 2*half_frame, 'EndPoints', 'discard');
    x = (half_frame:length(data)-half_frame)';
    % Linear trend in variance
    EQ = polyfit(x, var0, 1);
    x_mod = (1:length(data))';
    mod = polyval(EQ, x_mod);
    % Standardize variance-corrected data
    radar_Z(:,i) = data./sqrt(abs(mod));
end

% Define the vertical resolution of the core data
core_res = 0.02;

% Define the cutoff depth for radar traces and find index of crossover
% depth
cutoff = 25;
depth_bott = floor(min([min(radar.depth(end,:)) cutoff]));

% Trim radar traces to cutoff depth and interpolate data to vertical scale
% of the firn cores
radarZ_interp = zeros(depth_bott/core_res+1, size(radar.data_stack, 2));
for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:core_res:radar.depth(end,i));
    radarZ_i = interp1(radar.depth(:,i), radar_Z(:,i), depth_interp, 'pchip');
    radarZ_interp(:,i) = radarZ_i(1:size(radarZ_interp, 1));
end

% If manual layer picks are present, transform them to same depth and
% vertical scale as the interpolated radar data
if isfield(radar, 'man_layers')
    man_interp = zeros(size(radarZ_interp));
    for i = 1:size(radar.data_stack, 2)
%         depth_interp = (0:core_res:radar.depth(end,i))';
        layer_idx = logical(radar.man_layers(:,i));
        layer_num = radar.man_layers(layer_idx,i);
        man_depth_i = radar.depth(layer_idx,i);
        depth_idx = round(man_depth_i/core_res) + 1;
        man_interp(depth_idx(man_depth_i<=cutoff),i) = layer_num(man_depth_i<=cutoff);
    end
    radar.man_layers = man_interp;
end

% Assign structure output depth to interpolated depths
radar.depth = (0:core_res:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~20 m)
radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);

% Clear unnecessary variables
clearvars -except file cores Ndraw radar horz_res core_res

%% Iterative radon transforms

depth_interval = 4;
dist_interval = 250;
depth_sz = round(0.5*depth_interval/core_res);
dist_sz = round(0.5*dist_interval/horz_res);

s_matrix = nan(size(radar.data_smooth));
ii = dist_sz+1:3:size(radar.data_smooth,2)-dist_sz;
jj = depth_sz+1:5:size(radar.data_smooth,1)-depth_sz;

for i = 1:length(ii)
    for j = 1:length(jj)
        data = radar.data_smooth(jj(j)-depth_sz:jj(j)+depth_sz,ii(i)-dist_sz:ii(i)+dist_sz);
        theta = 0:179;
        [R,~] = radon(data, theta);
        
        R_max = sum(R.*(R>0));
        [~,theta_idx] = max(R_max);
        theta_max = theta(theta_idx);
%         [~,idx] = max(R, [], 2);
%         R_max(R_max<0) = 0;
%         weights = R_max/(sum(R_max));
%         theta_max = sum(weights.*theta);
%         theta_max = sum(weights.*theta(idx));
        s_matrix(jj(j),ii(i)) = -1*tand(theta_max-90);  
    end
    s_matrix(1,ii(i)) = 0;
    p = robustfit(1:length(radar.depth), s_matrix(:,ii(i)));  % Relation may be nonlinear
    s_matrix(end,ii(i)) = length(radar.depth)*p(2);
end


s_matrix(:,1) = s_matrix(:,ii(1));
s_matrix(:,end) = s_matrix(:,ii(end));


x = [1 ii size(radar.data_smooth,2)];
y = [1 jj size(radar.data_smooth,1)];
[X,Y] = meshgrid(radar.dist(x), radar.depth(y));
ss = s_matrix(y,x);
[Vx, Vy] = meshgrid(radar.dist, radar.depth);
Vq = interp2(X, Y, ss, Vx, Vy);
grad_smooth = imgaussfilt(Vq, 5);

% Diagnostic plot
ystart = 1:25:size(grad_smooth,1);
xstart = ones(1, length(ystart));
XY_raw = stream2(ones(size(grad_smooth)), grad_smooth, xstart, ystart, 1);
XY = XY_raw;
for k = 1:length(XY)
    XY{k}(:,1) = XY_raw{k}(:,1)*mean(diff(radar.dist));
    XY{k}(:,2) = XY_raw{k}(:,2)*.02;
end
figure
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
hold on
hlines = streamline(XY);
set(hlines, 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '--')
hold off


%% Find depth, width, and prominence of peaks for each radar trace

% Preallocate arrays for various components
peaks_raw = zeros(size(radar.data_smooth));
peak_width = zeros(size(radar.data_smooth));
Proms = cell(1, size(radar.data_smooth, 2));
widths = cell(1,size(radar.data_smooth, 2));
depths = cell(1, size(radar.data_smooth, 2));
depth_idx = cell(1, size(radar.data_smooth, 2));

for i = 1:size(radar.data_smooth, 2)
    data_i = radar.data_smooth(:,i);
    
    % Prominence threshold for peaks
    minProm = 0.50;
    
    % Min distance between peaks (in meters)
    minDist = 0.08;
    
    % Find peak statistics in each trace based on criteria
    [~, peaks_idx_i, widths_i, Prom_i] = findpeaks(data_i, ...
        'MinPeakProminence', minProm, ...
        'MinPeakDistance', minDist/core_res, 'WidthReference', 'halfheight');

    % Add peak prominence and width values to relevent matrices
    peaks_raw(peaks_idx_i,i) = Prom_i;
    peak_width(peaks_idx_i,i) = widths_i;
    
    % Add values to relevent cell arrays
    Proms{i} = Prom_i;
    widths{i} = widths_i;
    depths{i} = radar.depth(peaks_idx_i);
    depth_idx{i} = peaks_idx_i;
end


%%

[peak_group, layers2] = find_layers2(peaks_raw, peak_width, ...
    grad_smooth, core_res, horz_res);

% Preallocate arrays for the matrix indices of members of each layer
layers_idx2 = cell(1,length(layers2));
peaks2 = zeros(size(peaks_raw));

% For loop to coerce layers to have one row position for each trace
for i = 1:length(layers_idx2)
    
    % Find matrix indices of all members of ith layer
    layer_i = layers2{i};
    
    % Find row and col indices of members of ith layer
    [row, col] = ind2sub(size(radar.data_smooth), layer_i);
    mag = peaks_raw(layer_i);
    
    % Interpolate data to all column positions within the range of the
    % layer
    col_interp = min(col):max(col);
    
    % Interpolate row positions using a cubic smoothing spline
    row_interp = round(fnval(csaps(col, row), col_interp));
    row_interp(row_interp < 1) = 1;
    row_interp(row_interp > size(peaks2,1)) = size(peaks2,1);
    
    % Interpolate peak prominence magnitudes to all columns in range using
    % a cubic smoothing spline
    mag_interp = csaps(col, mag, 1/length(col_interp), col_interp);
    
    % Assign interpolated layer to output
    layer_interp = sub2ind(size(peaks2), row_interp, col_interp);
    peaks2(layer_interp) = mag_interp;
    layers_idx2{i} = layer_interp';
end

% Create matrix of layer group assignments
group_num = zeros(size(peaks2));
for i = 1:length(layers_idx2)
    group_num(layers_idx2{i}) = i;
end

%%

% Calculate continuous layer distances for each layer (accounting for 
% lateral size of stacked radar trace bins)
layers_dist = cellfun(@(x) numel(x)*horz_res, layers_idx2);

% Map layer prominence-distance values to the location within the radar
% matrix of the ith layer
layer_peaks = zeros(size(peaks2));
for i = 1:length(layers_idx2)
    layer_peaks(layers_idx2{i}) = peaks2(layers_idx2{i}).*layers_dist(i);
end


% Output layer arrays to radar structure
radar.peaks = peaks2;
radar.layers = layers_idx2;
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
    
    % Assign the 50% likelihood point based on median trace prominence and
    % layer length
    P_50 = median(Proms{i})*min([5000 0.5*radar.dist(end)]);
%     P_50 = median(Proms{i})*mean(cellfun(@length, layers_idx));
    
    % Assign min/max layer likelihoods, and calculate the logistic rate
    % coefficient
    Po = 0.05;
    K = 1;
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = layer_peaks(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar.depth(peaks_idx);
    
    % Likelihood of layer representing a year based on a logistic function
    % with rate (r) calculated above
    likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i));
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

end