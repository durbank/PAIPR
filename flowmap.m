% Script to generate the necessary figures used in the methods flowmap for
% AGU 2018


% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'work';
        %         computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
                
            case 'laptop'
                data_path = 'F:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/durba/Documents/MATLAB/Addons/';
        end
        
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = fullfile(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))
% Add export_fig to path
addon_folder = fullfile(addon_path, 'altmany-export_fig-cafc7c5/');
addpath(genpath(addon_folder))
% Add core functions to path
addpath(genpath('core'))

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);
Ndraw = 100;

% output_dir = uigetdir(data_path, ...
%     'Select directory to which to output images');

%%


tmp = radar_format(fullfile(data_path, ...
    '/radar/SEAT_Traverses/SEAT2010Kuband/flow-map'));
radar = tmp.segment(1);

% Convert to depth
[radar] = radar_depth(radar, cores);

% Find the mean response with depth in the radar data attributes across a
% given horizontal resolution (in meters)
horz_res = 25;
[radar] = radar_stack(radar, horz_res);

T = 710;
figure
imagesc(radar.dist, radar.depth(:,T), radar.data_stack)
hold on
plot([radar.dist(T) radar.dist(T)], [min(radar.depth(:,T)) max(radar.depth(:,T))], ...
    'b', 'LineWidth', 2)
xlabel('Distance (m)')
ylabel('Depth (m)')
ylim([0 25])
hold off





% Stationarize the radar response by differencing traces with a smoothing 
% spline
s = zeros(size(radar.data_stack));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth(:,i), radar.data_stack(:,i), 0.95, radar.depth(:,i));
end

figure
plot(radar.data_stack(:,T), radar.depth(:,T), 'k')
set(gca, 'Ydir', 'reverse')
hold on
plot(s(:,T), radar.depth(:,T), 'b', 'LineWidth', 2)
ylabel('Depth (m)')
ylim([0 24])
hold off





radar_stat = radar.data_stack - s;
depth_raw = radar.depth(:,T);

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

% Assign structure output depth to interpolated depths
radar.depth = (0:core_res:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~20 m)
radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);

figure
hold on
plot(zscore(radar_stat(:,T)), depth_raw, 'b--')
set(gca, 'Ydir', 'reverse')
h2 = plot(radar.data_smooth(:,T), radar.depth, 'c');
h2.Color(4) = 0.80;
xlabel('Z-statistic')
ylabel('Depth (m)')
ylim([0 25])
hold off

% Clear unnecessary variables
% clearvars -except file cores Ndraw radar horz_res core_res

h_sz = round((250/25)/2);
figure
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
hold on
plot([radar.dist(T-h_sz) radar.dist(T-h_sz)], [radar.depth(601) radar.depth(800)], 'r')
plot([radar.dist(T+h_sz) radar.dist(T+h_sz)], [radar.depth(601) radar.depth(800)], 'r')
plot([radar.dist(T-h_sz) radar.dist(T+h_sz)], [radar.depth(601) radar.depth(601)], 'r')
plot([radar.dist(T-h_sz) radar.dist(T+h_sz)], [radar.depth(800) radar.depth(800)], 'r')
scatter(radar.dist(T), radar.depth(700), 40, 'rx', 'LineWidth', 2.5)
ylim([0 25])
hold off




data = radar.data_smooth(601:800,T-h_sz:T+h_sz);
data_depth = radar.depth(601:800);
data_dist = 0:25:25*(size(data,2)-1);

% Angles over which to perform radon transform
theta = 0:179;

% Calculate transform
[R,~] = radon(data, theta);

% Find the sum of all positive transformed values for each angle
R_sum = sum(R.*(R>0));

% Find the maximum positive sum value, and assign it's
% corresponding angle as the dominant angle of layer gradients
[~,max_idx] = max(R_sum);
theta_max = theta(max_idx);

% Transform dominant angle to true layer gradient and assign to
% preallocated matrix
layer_grad = repmat(-1*tand(theta_max-90), length(data_depth), length(data_dist));

% Diagnostic plot
ystart = -100:15:length(data_depth);
xstart = ones(1, length(ystart));
XY_raw = stream2(ones(size(layer_grad)), layer_grad, xstart, ystart, 1);
XY = XY_raw;
for k = 1:length(XY)
    XY{k}(:,1) = XY_raw{k}(:,1)*horz_res;
    XY{k}(:,2) = XY_raw{k}(:,2)*core_res+data_depth(1);
end


figure
imagesc(data_dist, data_depth, data)
hold on
hlines = streamline(XY);
set(hlines, 'Color', 'r', 'LineStyle', '--')
scatter(data_dist(round(length(data_dist)/2)), data_depth(round(length(data_depth)/2)),...
    50, 'rx', 'LineWidth', 2.5)
hold off










% Define depth/distance intervals over which to perform radon transforms
% (in meters), and calculate data matrix window size (in data bins)
depth_interval = 4;
dist_interval = 250;
depth_sz = round(0.5*depth_interval/core_res);
dist_sz = round(0.5*dist_interval/horz_res);

% Define (overlapping) window center points for iterative radon transforms
ii = dist_sz+1:4:size(radar.data_smooth,2)-dist_sz;
jj = depth_sz+1:round(depth_sz/10):size(radar.data_smooth,1)-depth_sz;

% Preallocate matrix for estimated layer slopes
s_matrix = nan(size(radar.data_smooth));

% Iteratively calculate local layer gradients (data array units) for
% overlapping local windows within the radargram
for i = 1:length(ii)
    for j = 1:length(jj)
        
        % Define the local data window
        data_i = radar.data_smooth(jj(j)-depth_sz:jj(j)+depth_sz,ii(i)-dist_sz:ii(i)+dist_sz);
        
        % Angles over which to perform radon transform
        theta = 0:179;
        
        % Calculate transform
        [R,~] = radon(data_i, theta);
        
        % Find the sum of all positive transformed values for each angle
        R_sum = sum(R.*(R>0));
        
        % Find the maximum positive sum value, and assign it's
        % corresponding angle as the dominant angle of layer gradients
        [~,max_idx] = max(R_sum);
        theta_max = theta(max_idx);

        % Transform dominant angle to true layer gradient and assign to
        % preallocated matrix
        s_matrix(jj(j),ii(i)) = -1*tand(theta_max-90);  
    end
    
    % Define the surface layer to have a gradient of 0
    s_matrix(1,ii(i)) = 0;
    
    % Extrapolate the layer gradients at the base of the radargram based on
    % the robust linear regression of the slopes of overlying layers
    p = robustfit(1:length(radar.depth), s_matrix(:,ii(i)));
    s_matrix(end,ii(i)) = length(radar.depth)*p(2);
end

% Extrapolate gradient values for radargram edges based on nearest non-NaN
% values
s_matrix(:,1) = s_matrix(:,ii(1));
s_matrix(:,end) = s_matrix(:,ii(end));

% Update the along-trace and depth locations of non-NaN values in slope
% matrix
x = [1 ii size(radar.data_smooth,2)];
y = [1 jj size(radar.data_smooth,1)];

% Interpolate layer gradients for every trace and depth point
[X,Y] = meshgrid(radar.dist(x), radar.depth(y));
ss = s_matrix(y,x);
[Vx, Vy] = meshgrid(radar.dist, radar.depth);
Vq = interp2(X, Y, ss, Vx, Vy);
grad_smooth = Vq;

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
ylim([0 25])
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


figure
imagesc(radar.dist, radar.depth, peaks_raw, [0 3.5])
hold on
hlines = streamline(XY);
set(hlines, 'Color', 'r', 'LineStyle', '--')
ylim([0 25])
hold off







%%

data_peaks = peaks_raw(601:800,T-h_sz:T+h_sz);
data_width = peak_width(601:800,T-h_sz:T+h_sz);











