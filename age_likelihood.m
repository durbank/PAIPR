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
    filesep 'ProcessedSEAT2010' filesep 'grid_SEAT10_3' filesep]);

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
P_50 = 750;
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

% Caculate radon-transform weighting coefficient for radar file based on
% average segment length
% w_RT = mean(layer_ima(layer_ima>0));
w_RT = K*Po./(Po + (K-Po)*exp(-r*mean(lengths)));

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


%% Calculate RT and peak likelihoods separately
age_peaks = zeros([size(radar.data_smooth) Ndraw]);
age_RT = zeros([size(radar.data_smooth) Ndraw]);
age_w = zeros([size(radar.data_smooth) Ndraw]);
data = radar.data_smooth;
depth_i = radar.depth_interp;
parfor i = 1:size(radar.data_smooth, 2)
    
    %%Calculate RT likelihoods
    likelihood_RT = layer_ima(:,i);
    
    %%Calculate peak likelihoods
    data_i = data(:,i);
    likelihood_peaks = zeros(length(data_i), 1);
    % minProm = 0.01;
    minProm = 0.1*iqr(data_i);
    [~, idx_peaks, ~, Prom] = findpeaks(data_i, 'MinPeakProminence', minProm);
    
    % Solve for r so that prominence with 0.50 probability is equal to
    % interquartile range of the smoothed data (may need to be optimized to
    % observations using a scaling constant)
    Po = 0.01;
    K = 1;
    P_50 = 1*iqr(data_i);
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % % Calculate distance weights based off of the distances to the nearest
    % % four peaks (in future should scale the effect based off the variance
    % % in the expected distance, dist_Ex_i)
    % [~, Prom_dist] = knnsearch(depths_peaks, depths_peaks, 'K', 5);
    % dist_Ex_i = f(1)*depths_peaks + f(2);
    % w = sum([0.5 0.5 0.25 0.25].*(Prom_dist(:,2:end)-dist_Ex_i), 2)./dist_Ex_i;
    % Pw = P.^(1-w);
    % P = Pw;
    
    % Probability of peak representing a year based on a logistic function
    % with rate (r) calculated above
    P_peaks = K*Po./(Po + (K-Po)*exp(-r*Prom));
    likelihood_peaks(idx_peaks) = P_peaks;
    
    %%%Compare age-depth profiles from each method and the weighted average
    
    %%Calculate age-depth profile distribution from peak likelihoods
    
    depths_peaks = depth_i(idx_peaks);
    yr_idx = zeros(length(depths_peaks), Ndraw);
    for j = 1:length(depths_peaks)
        R = rand(Ndraw, 1) <= P_peaks(j);
        yr_idx(j,:) = R;
    end
    age_i = zeros(length(data_i), Ndraw);
    for j = 1:Ndraw
        depths_j = [0; depths_peaks(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        %     age_peaks(:,j) = interp1(depths_j, yrs_j, depth_i, 'linear', 'extrap');
        age_i(:,j) = interp1(depths_j, yrs_j, depth_i, 'linear', 'extrap');
    end
    age_peaks(:,i,:) = age_i;
    
    %%Calculate age-depth profile distribution from RT likelihoods
    
    depths_RT = depth_i(layer_ima(:,i)>0);
    P_RT = layer_ima(:,i);
    P_RT = P_RT(P_RT>0);
    yr_idx = zeros(length(depths_RT), Ndraw);
    for j = 1:length(depths_RT)
        R = rand(Ndraw, 1) <= P_RT(j);
        yr_idx(j,:) = R;
    end
    age_i = zeros(length(data_i), Ndraw);
    for j = 1:Ndraw
        depths_j = [0; depths_RT(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        age_i(:,j) = interp1(depths_j, yrs_j, depth_i, 'linear', 'extrap');
    end
    age_RT(:,i,:) = age_i;
    
    %%Shift RT likelihood depths to match nearby depths in peak likelihood
    
    tol = 3;
    idx_RT = find(layer_ima(:,i));
    RT_shift = zeros(length(likelihood_RT), 1);
    for j = idx_RT'
        if j<=tol
            data_j = likelihood_peaks(1:j+tol);
        elseif j+tol>length(likelihood_peaks)
            data_j = likelihood_peaks(j-tol:end);
        else
            data_j = likelihood_peaks(j-tol:j+tol);
        end
        if isempty(data_j(data_j>0))
            RT_shift(j) = likelihood_RT(j);
        else
            data_j(data_j<=0) = -1;
            [~, min_idx] = min(abs(data_j - likelihood_RT(j)));
            RT_shift(j-tol+min_idx-1) = likelihood_RT(j);
        end
    end
    
    %%Calculate age-depth profile distribution from weighted likelihoods
    
    % Combine two likelihoods into a single weighted average
    likelihood_w = (likelihood_peaks + w_RT*RT_shift)/(1+w_RT);
    
    depths_weighted = depth_i(likelihood_w>0);
    P_w = likelihood_w(likelihood_w>0);
    yr_idx = zeros(length(depths_weighted), Ndraw);
    for j = 1:length(depths_weighted)
        R = rand(Ndraw, 1) <= P_w(j);
        yr_idx(j,:) = R;
    end
    age_i = zeros(length(data_i), Ndraw);
    for j = 1:Ndraw
        depths_j = [0; depths_weighted(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        age_i(:,j) = interp1(depths_j, yrs_j, depth_i, 'linear', 'extrap');
    end
    age_w(:,i,:) = age_i;
end

% Diagnostic figure
figure
hold on
h1 = plot(radar.depth_interp, mean(mean(age_peaks, 3), 2), 'm');
h2 = plot(radar.depth_interp, mean(mean(age_RT, 3), 2), 'c');
h3 = plot(radar.depth_interp, mean(mean(age_w, 3), 2), 'r');
plot(radar.depth_interp, mean(mean(age_w, 3), 2)+std(mean(age_w, 3), [], 2), 'r--')
plot(radar.depth_interp, mean(mean(age_w, 3), 2)-std(mean(age_w, 3), [], 2), 'r--')
h4 = plot(core.depth, core.age, 'b');
legend([h1 h2 h3 h4], 'Mean peak age', 'Mean RT age', 'Mean weighted age', 'Core age')
hold off

% % Diagnostic figure
% figure
% hold on
% plot(depth_i, mean(age_peaks, 2), 'b')
% plot(depth_i, mean(age_peaks, 2) + std(age_peaks, [], 2), 'b--')
% plot(depth_i, mean(age_peaks, 2) - std(age_peaks, [], 2), 'b--')
% plot(depth_i, mean(age_RT, 2), 'r')
% plot(depth_i, mean(age_RT, 2) + std(age_RT, [], 2), 'r--')
% plot(depth_i, mean(age_RT, 2) - std(age_RT, [], 2), 'r--')
% plot(depth_i, mean(age_w, 2), 'm', 'LineWidth', 2)
% plot(depth_i, mean(age_w, 2) + std(age_w, [], 2), 'm--')
% plot(depth_i, mean(age_w, 2) - std(age_w, [], 2), 'm--')
% plot(core.depth, core.age, 'k', 'LineWidth', 2)
% hold off
