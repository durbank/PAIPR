% This is an alpha script to combine the peak likelihood estimates and 
% radon transform estimates to produce a single age-depth profile 
% distribution

% Assigns the path to the data and add-on folders, based on whether running
% from Durban's Macbook Pro or a lab PC
PC_true = ispc;
switch PC_true
    case true
        data_path = 'D:/Research/Antarctica/WAIS Variability/';
        addon_path = 'D:/Research/Antarctica/WAIS Variability/Addons/';
    case false
        data_path = '/Volumes/WARP/Research/Antarctica/WAIS Variability/';
        addon_path = '/Users/Durbank/Documents/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

% Import firn core data
[cores] = import_cores(strcat(data_path, ['SEAT_cores' filesep ...
    'DGK_core_data.xlsx']));

%%
% Directory to radar files of interest
radar_dir = strcat(data_path, ['SEAT_Traverses' filesep 'SEAT2010Kuband'...
    filesep 'ProcessedSEAT2010' filesep 'grid_SEAT10_6' filesep]);

% List all files matching 'wild' within radar directory
wild = 'layers*';
files = dir(strcat(radar_dir, wild));

% Select individual radar file for analysis
i = randi(length(files));
file = strcat([files(i).folder filesep], files(i).name);

% Number of simulations to perform on age-depth Monte Carlo
Ndraw = 100;

% Conversion to depth
[radar, core] = radar_depth(file, cores);

% Stationarize the radar response using a smoothing spline
s = zeros(size(radar.data_out));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth, radar.data_out(:,i), 0.95, radar.depth);
end
radar_stat = radar.data_out - s;

% Remove linear trend in variance (attentuation with depth) and convert to
% z-score statistics
radar_Z = zeros(size(radar_stat));
for i = 1:size(radar_stat, 2)
    data = radar_stat(:,i);
    half_frame = round(0.5*length(data)/5);
    var0 = movvar(data, 2*half_frame, 'EndPoints', 'discard');
    x = (half_frame:length(data)-half_frame)';
    EQ = polyfit(x, var0, 1);
    x_mod = (1:length(data))';
    mod = EQ(1)*x_mod+EQ(2);
    radar_Z(:,i) = zscore(data./sqrt(abs(mod)));
end

% Resample and interpolate radar depth and data to a constant depth
% resolution using a shape-preserving cubic interpolation
resolution = 0.02;       % Depth resolution in meters
radar.depth_interp = (0:resolution:radar.depth(end))';
radarZ_resamp = interp1(radar.depth, radar_Z, radar.depth_interp, 'pchip');

% Find the mean response with depth in the resampled radar data across a
% given lateral distance 'window' (lateral bin size, where 1 bin ~ 2.5 m)
window = 50;
radar_mean = movmean(radarZ_resamp, window, 2);

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~20 m)
radar.data_smooth = sgolayfilt(radar_mean, 3, 9);

%%
%Add radon scripts to path
addpath LTARE_codes/

% Generate binary image of annual layers
[CC, data_pts] = DGK_horizons(radar.data_smooth);

% Generate matrix of annual horizons
L = labelmatrix(CC);

% Statistics on each horizon segment
stats = regionprops(CC, 'PixelList');

% Define variables for the likelihood logistic function
Po = 0.01;
K = 1;
P_50 = 500;
r = log((K*Po/0.50-Po)/(K-Po))/-P_50;

% Preallocate weight matrix and stats structure
layer_ima = zeros(size(radar.data_smooth));
stats_new = struct('PixelList', cell(1, length(stats)), ...
    'Layer_length', [], 'Likelihood', []);

for i = 1:length(stats)
    % Indices of unique x values within radar data matrix for ith layer
    x_all = unique(stats(i).PixelList(:,1));
    
    % Boundaries between layer segments with the same x values
    int_jump = [1; diff(stats(i).PixelList(:,1))];
    x_idx = [find(int_jump); length(stats(i).PixelList(:,1))];
    
    % Find the length of the continuous layer (in meters)
    layer_length = radar.dist(x_all(end)) - radar.dist(x_all(1));
    stats_new(i).Layer_length = layer_length;
    
    % Assign layer likelihood based on logistic function
    likelihood = K*Po./(Po + (K-Po)*exp(-r*layer_length));
    stats_new(i).Likelihood = likelihood;
    
    y_all = zeros(size(x_all));
    for j = 1:length(x_all)
        % Average all y-values for each unique x in ith layer (coerces the
        % layer to be 1D)
        y_all(j) = round(mean(stats(i).PixelList(x_idx(j):x_idx(j+1),2)));
        
        % Add layer likelihoods to layer image matrix
        layer_ima(y_all(j),x_all(j)) = likelihood;
    end
    
    % Add 1D layer indices to stats structure
    stats_new(i).PixelList = [x_all y_all];
end

% Vector of layer segment lengths
lengths = extractfield(stats_new, 'Layer_length');

% Vector of layer segment likelihoods
P = extractfield(stats_new, 'Likelihood');

% Plot segment lengths vs. segment likelihoods
figure
plot(lengths, P, 'o')

% % Caculate radon-transform weighting coefficient for radar file based on 
% % average segment length
w_RT = mean(layer_ima(layer_ima>0));

%%
% Year associated with the first pick of the algorithm
% age_top = round(radar.collect_date);
age_top = radar.collect_date;
yr_pick1 = ceil(age_top-1);

% Indices of annual horizions in core
yr_idx = logical([diff(floor(core.age)); 0]);

% Difference in depth between annual horizions in core
depth_diff = diff(core.depth(yr_idx));
depth_btw_yr = [depth_diff; mean(depth_diff)];

% Model of distance between annual horizions with depth in core
f = polyfit(core.depth(yr_idx), depth_btw_yr, 1);

% Expected thickness (and standard deviation) of annual layers 
% throughout the depth of the core
dist_Ex = f(1)*radar.depth_interp + f(2);
dist_STD = std(depth_btw_yr);


%% PURE TESTING
i = 100;
likelihood_peaks = zeros(size(radar.data_smooth, 1), 1);

% Estimate the depth of annual peaks within each smoothed radar trace
data_i = radar.data_smooth(:,i);
depth_i = radar.depth_interp;
minProm = 0.01;
minProm = 0.1*iqr(data_i);
[~, idx_peaks, ~, Prom] = findpeaks(data_i, 'MinPeakProminence', minProm);

% Solve for r so that prominence with 0.50 probability is equal to
% interquartile range of the smoothed data (may need to be optimized to
% observations using a scaling constant)
Po = 0.01;
K = 1;
P_50 = 1*iqr(data_i);
r = log((K*Po/0.50-Po)/(K-Po))/-P_50;

% Probability of peak representing a year based on a logistic function
% with rate (r) calculated above
P_peaks = K*Po./(Po + (K-Po)*exp(-r*Prom));
likelihood_peaks(idx_peaks) = P_peaks;

% RT likelihoods
%depth_RT = depth_i(layer_ima(:,i)>0);
likelihood_RT = layer_ima(:,i);

% Diagnostic figure
figure
plot(depth_i, likelihood_peaks, 'LineWidth', 2)
hold on
plot(depth_i, likelihood_RT)
hold off

%%Compare age-depth profiles from the two independent methods
depths_peaks = depth_i(idx_peaks);

yr_idx = zeros(length(depths_peaks), Ndraw);
for j = 1:length(depths_peaks)
    R = rand(Ndraw, 1) <= P_peaks(j);
    yr_idx(j,:) = R;
end

age_peaks = zeros(length(data_i), Ndraw);
for j = 1:Ndraw
    depths_j = [0; depths_peaks(logical(yr_idx(:,j)))];
    yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
    age_peaks(:,j) = interp1(depths_j, yrs_j, depth_i, 'linear', 'extrap');
end


depths_RT = depth_i(layer_ima(:,i)>0);
P_RT = layer_ima(:,i);
P_RT = P_RT(P_RT>0);

yr_idx = zeros(length(depths_RT), Ndraw);
for j = 1:length(depths_RT)
    R = rand(Ndraw, 1) <= P_RT(j);
    yr_idx(j,:) = R;
end

age_RT = zeros(length(data_i), Ndraw);
for j = 1:Ndraw
    depths_j = [0; depths_RT(logical(yr_idx(:,j)))];
    yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
    age_RT(:,j) = interp1(depths_j, yrs_j, depth_i, 'linear', 'extrap');
end

figure
hold on
plot(depth_i, mean(age_peaks, 2))
plot(depth_i, mean(age_RT, 2))
hold off
    




%%


% Estimate the depth of annual peaks within each smoothed radar trace
ages = zeros([size(radar.data_smooth) Ndraw]);
data = radar.data_smooth;
depth_i = radar.depth_interp;
parfor i = 1:size(radar.data_smooth, 2)
    data_i = data(:,i);
    minProm = 0.01;
    minDist = abs(min(dist_Ex) - 3*dist_STD);
    minDist = 0;
    [~, depths_peaks, ~, Prom] = findpeaks(data_i, depth_i, ...
        'MinPeakProminence', minProm, 'MinPeakDistance', minDist);
    
    % Remove the surface horizon peak from the data
    [~,max_idx] = max(Prom);
    Prom(max_idx) = [];
    depths_peaks(max_idx) = [];
    
    % Solve for r so that prominence with 0.50 probability is equal to 
    % interquartile range of the smoothed data (may need to be optimized to
    % observations using a scaling constant)
    Po = 0.01;
    K = 1;
    P_50 = 1*iqr(data_i);
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % Probability of peak representing a year based on a logistic function
    % with rate (r) calculated above
    P = K*Po./(Po + (K-Po)*exp(-r*Prom));
    
%     % Function using a more general logistic function (requires
%     % optimization of k)
%     L = 1;
%     P_50 = mean(Prom);
%     k = 4;
%     P = L./(1 + exp(-k.*(Prom-P_50)));

    % Calculate distance weights based off of the distances to the nearest
    % four peaks (in future should scale the effect based off the variance 
    % in the expected distance, dist_Ex_i)
    [~, Prom_dist] = knnsearch(depths_peaks, depths_peaks, 'K', 5);
    dist_Ex_i = f(1)*depths_peaks + f(2);
    w = sum([0.5 0.5 0.25 0.25].*(Prom_dist(:,2:end)-dist_Ex_i), 2)./dist_Ex_i;
    Pw = P.^(1-w);
    P = Pw;
    
    yr_idx = zeros(length(depths_peaks), Ndraw);
    for j = 1:length(depths_peaks)
        R = rand(Ndraw, 1) <= P(j);
        yr_idx(j,:) = R;
    end
    
    age_i = zeros(length(data_i), Ndraw);
    for j = 1:Ndraw
        depths_j = [0; depths_peaks(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        age_i(:,j) = interp1(depths_j, yrs_j, depth_i, 'linear', 'extrap');
    end
    ages(:,i,:) = age_i;
end
radar.age = ages;

radar = rmfield(radar, {'time_trace', 'arr_layers', 'TWTT'});


