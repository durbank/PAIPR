% Script to generate the necessary figures used in the methods flowmap for
% AGU 2018


% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'laptop';
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

output_dir = uigetdir(data_path, ...
    'Select directory in which to output images');

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
fig1 = figure;
imagesc(radar.dist, radar.depth(:,T), radar.data_stack)
hold on
plot([radar.dist(T) radar.dist(T)], [min(radar.depth(:,T)) max(radar.depth(:,T))], ...
    'b', 'LineWidth', 2)
xlabel('Distance (m)')
ylabel('Depth (m)')
ylim([0 25])
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [12, 7])
hold off

f1_nm = 'Radargram_raw';
export_fig(fig1, fullfile(output_dir, f1_nm), '-pdf', '-q101', '-cmyk');
close(fig1)



% Stationarize the radar response by differencing traces with a smoothing 
% spline
s = zeros(size(radar.data_stack));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth(:,i), radar.data_stack(:,i), 0.95, radar.depth(:,i));
end

fig2 = figure;
plot(radar.data_stack(:,T), radar.depth(:,T), 'k')
set(gca, 'Ydir', 'reverse')
hold on
plot(s(:,T), radar.depth(:,T), 'b', 'LineWidth', 2)
ylabel('Depth (m)')
ylim([0 24])
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 7])
hold off

f2_nm = 'stationarize';
export_fig(fig2, fullfile(output_dir, f2_nm), '-pdf', '-q101', '-cmyk');
close(fig2)




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

fig3 = figure;
hold on
plot(zscore(radar_stat(:,T)), depth_raw, 'c--')
set(gca, 'Ydir', 'reverse')
h2 = plot(radar.data_smooth(:,T), radar.depth, 'b');
h2.Color(4) = 0.80;
xlabel('Z-statistic')
ylabel('Depth (m)')
xlim([-3 4])
ylim([0 25])
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 7])
hold off

f3_nm = 'trace_smooth';
export_fig(fig3, fullfile(output_dir, f3_nm), '-pdf', '-q101', '-cmyk');
close(fig3)

h_sz = round((250/25)/2);
data = radar.data_smooth(601:800,T-h_sz:T+h_sz);
data_depth = radar.depth(601:800);
data_dist = 0:25:25*(size(data,2)-1);
row_n = 122;
col_n = 6;




fig4 = figure;
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
hold on
plot([radar.dist(T-h_sz) radar.dist(T-h_sz)], [radar.depth(601) radar.depth(800)],...
    'r', 'LineWidth', 2)
plot([radar.dist(T+h_sz) radar.dist(T+h_sz)], [radar.depth(601) radar.depth(800)],...
    'r', 'LineWidth', 2)
plot([radar.dist(T-h_sz) radar.dist(T+h_sz)], [radar.depth(601) radar.depth(601)],...
    'r', 'LineWidth', 2)
plot([radar.dist(T-h_sz) radar.dist(T+h_sz)], [radar.depth(800) radar.depth(800)],...
    'r', 'LineWidth', 2)
scatter(radar.dist(T-h_sz+col_n-1), radar.depth(601+row_n), 50, 'r*')
xlabel('Distance (m)')
ylabel('Depth (m)')
ylim([0 25])
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [12, 7])
hold off

f4_nm = 'radargram_smooth';
export_fig(fig4, fullfile(output_dir, f4_nm), '-pdf', '-q101', '-cmyk');
close(fig4)


%%

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


fig5 = figure;
imagesc(data_dist, data_depth, data)
hold on
hlines = streamline(XY);
set(hlines, 'Color', 'r', 'LineStyle', '--')
scatter(data_dist(col_n), data_depth(row_n), 50, 'r*')
xlabel('Distance (m)')
ylabel('Depth (m)')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 7])
hold off

fig_nm = 'RT_subset';
export_fig(fig5, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk');
close(fig5)









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

fig6 = figure;
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
hold on
hlines = streamline(XY);
set(hlines, 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '--')
plot([radar.dist(T-h_sz) radar.dist(T-h_sz)], [radar.depth(601) radar.depth(800)], ...
    'm', 'LineWidth', 2)
plot([radar.dist(T+h_sz) radar.dist(T+h_sz)], [radar.depth(601) radar.depth(800)],...
    'm', 'LineWidth', 2)
plot([radar.dist(T-h_sz) radar.dist(T+h_sz)], [radar.depth(601) radar.depth(601)],...
    'm', 'LineWidth', 2)
plot([radar.dist(T-h_sz) radar.dist(T+h_sz)], [radar.depth(800) radar.depth(800)],...
    'm', 'LineWidth', 2)
scatter(radar.dist(T-h_sz+col_n-1), radar.depth(601+row_n), 50, 'm*')
xlabel('Distance (m)')
ylabel('Depth (m)')
ylim([0 25])
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [12, 7])
hold off

fig_nm = 'radargram_RT';
export_fig(fig6, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk');
close(fig6)





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


% figure('Position', [50 50 1500 800])
% imagesc(radar.dist, radar.depth, peaks_raw, [0 3.5])
% hold on
% hlines = streamline(XY);
% set(hlines, 'Color', 'r', 'LineStyle', '--')
% xlabel('Distance (m)')
% ylabel('Depth (m)')
% ylim([0 25])
% hold off







%%

data_peaks = peaks_raw(601:800,T-h_sz:T+h_sz);
data_width = peak_width(601:800,T-h_sz:T+h_sz);

row_n = 122;
col_n = 6;
n_idx = sub2ind(size(data_peaks), row_n, col_n);
mag_i = data_peaks(n_idx);
width_i = data_width(n_idx);

peak_local = data_peaks;
peak_local(n_idx) = 0;
local_idx = find(peak_local);
mag_local = peak_local(local_idx);
[row_local, col_local] = ind2sub(size(peak_local), local_idx);
data_local = [row_local col_local];
% data_local = [(row_idx(1)+row_local-1) (col_idx(1)+col_local-1)];

stream_n1 = stream2(ones(size(data_peaks)), layer_grad, ...
    col_n, row_n, [0.10, 1000]);
stream_n2 = stream2(-ones(size(data_peaks)), -layer_grad, ...
    col_n, row_n, [0.10, 1000]);
cols_stream = vertcat(flipud(stream_n2{1}(:,1)), stream_n1{1}(:,1));
rows_stream = vertcat(flipud(stream_n2{1}(:,2)), stream_n1{1}(:,2));
[~,u_idx] = unique(cols_stream);
vq = interp1(cols_stream(u_idx), rows_stream(u_idx), 1:size(data_peaks, 2), ...
    'linear', 'extrap');
data_stream = [vq' (1:size(data_peaks,2))'];


dist_n = sqrt(((data_local(:,1)-(data_stream(:,1))')./...
            (0.5*width_i)).^2 + (data_local(:,2)-(data_stream(:,2))').^2 + ...
            (repmat(mag_local-mag_i, 1, size(data_stream,1))).^2);
threshold = 3;
dist_idx = min(dist_n, [], 2) <= threshold;


fig7 = figure;
imagesc(data_dist, data_depth, data_peaks)
hold on
plot(horz_res*(cols_stream-1), data_depth(1)+(rows_stream-1)*core_res, 'r--')
scatter(data_dist(col_n), data_depth(row_n), 50, 'm*')
scatter(data_dist(col_local(dist_idx)), data_depth(row_local(dist_idx)),...
    25, 'filled', 'm')
scatter(data_dist(col_local(~dist_idx)), data_depth(row_local(~dist_idx)), 25, 'bx')
xlabel('Distance (m)')
ylabel('Depth (m)')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 7])
hold off

fig_nm = 'neighbors';
export_fig(fig7, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk');
close(fig7)



%%

% Find continuous layers within radargram based on peaks and layer stream
% field
[peak_group, layers] = find_layers2(peaks_raw, peak_width, ...
    grad_smooth, core_res, horz_res);

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



fig8 = figure;
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
hold on
for i = 1:length(radar.layers)
    [row,col] = ind2sub(size(radar.data_smooth), radar.layers{i});
    plot(radar.dist(col), radar.depth(row), 'LineWidth', 2)
end
plot([radar.dist(T-h_sz) radar.dist(T-h_sz)], [radar.depth(601) radar.depth(800)], ...
    'Color', [0.90 0.40 0.15], 'LineWidth', 2.5)
plot([radar.dist(T+h_sz) radar.dist(T+h_sz)], [radar.depth(601) radar.depth(800)],...
    'Color', [0.90 0.40 0.15], 'LineWidth', 2.5)
plot([radar.dist(T-h_sz) radar.dist(T+h_sz)], [radar.depth(601) radar.depth(601)],...
    'Color', [0.90 0.40 0.15], 'LineWidth', 2.5)
plot([radar.dist(T-h_sz) radar.dist(T+h_sz)], [radar.depth(800) radar.depth(800)],...
    'Color', [0.90 0.40 0.15], 'LineWidth', 2.5)
ylim([0 25])
xlabel('Distance (m)')
ylabel('Depth (m)')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [12, 7])
hold off

fig_nm = 'radar_layers';
export_fig(fig8, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk');
close(fig8)


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
    
%     % Assign the 50% likelihood point based on median trace prominence and
%     % layer length
%     P_50 = median(Proms{i})*max([10000 0.25*radar.dist(end)]);
%     
%     % Assign min/max layer likelihoods, and calculate the logistic rate
%     % coefficient
%     Po = 0.05;
%     K = 1;
%     r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = layer_peaks(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar.depth(peaks_idx);
    
    % Likelihood of layer representing a year based on a logistic function
    % with rate (r) calculated above
    r = -2.4333e-4; % [-3.18e-4 -1.55e-4]
    k = 4.4323;     % [3.25 4.8]
    
    likelihood = 1./(1+exp(r*peaks_i + k));
%     likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i));
    radar.likelihood(peaks_idx,i) = likelihood;
    
    % Assign MC simulation annual layer presence based on layer likelihood
    % values
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end

%     yr_idx = zeros(length(depths_i), Ndraw);
%     for j = 1:length(depths_i)
%         yr_idx(j,:) = likelihood(j) >= 0.5;
%     end    

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


%%

data_intdist = layer_peaks(601:800,T-h_sz:T+h_sz);
dist_peaks = data_intdist(data_intdist>0);
likelihood = 1./(1+exp(r*dist_peaks + k));

fig9 = figure;
imagesc(data_dist, data_depth, data_intdist)
xlabel('Distance (m)')
ylabel('Depth (m)')
xlabel('Distance (m)')
ylabel('Depth (m)')
c = colorbar;
c.Label.String = 'Integrated prominence-distance';
c.Label.FontSize = 12;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 7])

fig_nm = 'layer_vals';
export_fig(fig9, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk');
close(fig9)


x = 0:max(dist_peaks);
y = 1./(1+exp(r*x + k));

fig10 = figure;
hold on
plot(x,y, 'k--')
scatter(dist_peaks, likelihood, 10, [0.90 0.40 0.15], 'filled')
xlabel('Integrated prominence-distance')
ylabel('Likelihood of annual layer')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [6, 7])
hold off

fig_nm = 'logistic';
export_fig(fig10, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk');
close(fig10)


fig11 = figure;
imagesc(radar.dist, radar.depth, mean(radar.ages,3))
hold on
plot([radar.dist(T) radar.dist(T)], [radar.depth(1) radar.depth(end)], ...
    'b', 'LineWidth', 2)
colormap('hsv')
c = colorbar;
c.Label.String = 'Calendar year';
c.Label.FontSize = 12;
caxis([1970 2010])
ylim([0 25])
xlabel('Distance (m)')
ylabel('Depth (m)')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [12, 7])
hold off

fig_nm = 'radar_ages';
export_fig(fig11, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk');
close(fig11)


fig12 = figure;
hold on
% for n = 1:100
%     h0 = plot(radar.ages(:,T,n), radar.depth, 'b', 'LineWidth', 0.5);
%     h0.Color(4) = 0.10;
% end
plot(mean(squeeze(radar.ages(:,T,:)),2), radar.depth, 'b', 'LineWidth', 2)
plot(mean(squeeze(radar.ages(:,T,:)),2) + 2*std(squeeze(radar.ages(:,T,:)),[],2), ...
    radar.depth, 'b--')
plot(mean(squeeze(radar.ages(:,T,:)),2) - 2*std(squeeze(radar.ages(:,T,:)),[],2), ...
    radar.depth, 'b--')
set(gca, 'Ydir', 'reverse')
xlim([1970 2012])
ylim([0 25])
xlabel('Distance (m)')
ylabel('Depth (m)')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 7])
hold off

fig_nm = 'trace_age';
export_fig(fig12, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk');
close(fig12)


