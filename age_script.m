% Directory to radar files of interest

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

% Calculate radar ages and associated other data
[radar, core] = age_MC(file, cores, Ndraw);

%%
% Automated core age
yr_top = core.age(1);
yr_pick1 = floor(yr_top);

iso = sgolayfilt((zscore(core.dD)+zscore(core.d18O))/2, 3, 9);
minDist = 0.02;    % At 2 cm resolution, = 02 cm min distance
minProm = 0.01;
[~, depths_peaks, ~, Prom] = findpeaks(iso, core.depth, ...
    'MinPeakDistance', minDist, 'MinPeakProminence', minProm);
Po = 0.01;
K = 1;
P_50 = iqr(iso);
r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
P = K*Po./(Po + (K-Po)*exp(-r*Prom));
% plot(Prom, P, 'ro')

yr_idx = zeros(length(depths_peaks), Ndraw);
for i = 1:length(depths_peaks)
    R = rand(Ndraw, 1) <= P(i);
    yr_idx(i,:) = R;
end

core.age_auto = zeros(length(core.depth), Ndraw);
for i = 1:Ndraw
    depths_i = [0; depths_peaks(logical(yr_idx(:,i)))];
    yrs_i = ([yr_top yr_pick1:-1:yr_pick1-length(depths_i)+2])';
    core.age_auto(:,i) = interp1(depths_i, yrs_i, core.depth, 'linear', 'extrap');
end

%%

% Select random radar trace for comparison plots
i = randi(size(radar.data_smooth, 2));

% Plot radargram
figure('Position', [200 200 1500 800])
imagesc([0 radar.dist(end)], [0 radar.depth_interp(end)], radar.data_smooth, [-2 2])
colorbar
xlabel('Distance along profile (m)')
ylabel('Depth (m)')
hold on
plot([radar.dist(i) radar.dist(i)], [0 radar.depth_interp(end)], 'r', 'LineWidth', 2)
xlim([0 radar.dist(end)])
ylim([0 radar.depth_interp(end)])
set(gca, 'Ydir', 'reverse', 'FontSize', 18)
hold off

% Compare radar signal at random trace to annual horizion selection at that
% trace
% layer_idx = logical([diff(floor(squeeze(radar.age(:,i,100)))); 0]);
layer_idx = logical([diff(floor(squeeze(radar.age(:,i,randi(Ndraw))))); 0]);
figure('Position', [50 50 500 1200])
hold on
plot(radar.data_smooth(:,i), radar.depth_interp, 'r', 'LineWidth', 1.5)
scatter(radar.data_smooth(layer_idx,i), radar.depth_interp(layer_idx), 'b<', 'filled')
xlim([-1.25 2])
xlabel('Radar Z-score')
ylabel('Depth (m)')
set(gca, 'Ydir', 'reverse')
hold off

% Compare radar signal at random trace to the mean annual horizions at that
% same random trace
layer_idx = logical([diff(floor(mean(squeeze(radar.age(:,i,:)), 2))); 0]);
figure('Position', [50 50 1200 500])
hold on
plot(radar.depth_interp, radar.data_smooth(:,i))
plot([radar.depth_interp(layer_idx) radar.depth_interp(layer_idx)], ...
    [-1.5 1.5], 'k');
ylim([-1.5 1.5])

% Compare core age to radar age
n_plot = floor(sqrt(Ndraw/3));
plot_idx = reshape(randi(size(radar.data_smooth, 2), n_plot), 1, n_plot^2);

figure
hold on
for j = 1:length(plot_idx)
    plot(radar.depth_interp, radar.age(:,plot_idx(j)), 'LineWidth', 0.25)
end
h1 = plot(core.depth, core.age, 'k', 'LineWidth', 3);
plot(core.depth, mean(core.age_auto, 2)+std(core.age_auto, 1, 2), 'r--', 'LineWidth', 2)
plot(core.depth, mean(core.age_auto, 2)-std(core.age_auto, 1, 2), 'r--', 'LineWidth', 2)
h2 = plot(core.depth, mean(core.age_auto, 2), 'r', 'LineWidth', 3);
xlabel('Depth (m)')
ylabel('Calendar Year')
legend([h1], 'Manual core picks')
hold off


plot_idx = randi(Ndraw, 1, round(Ndraw/3));
figure
hold on
for j = 1:length(plot_idx)
    plot(radar.depth_interp, radar.age(:,i,plot_idx(j)), 'LineWidth', 0.25)
end
h1 = plot(core.depth, core.age, 'k', 'LineWidth', 3);
h2 = plot(radar.depth_interp, mean(squeeze(radar.age(:,i,:)), 2), 'r', 'LineWidth', 3);
xlabel('Depth (m)')
ylabel('Calendar Year')
legend([h1 h2], 'Manual core picks', 'Mean radar age')
hold off


% Plot core vs mean radar ages (and uncertainty) for whole file
age_mean = mean(mean(radar.age, 3), 2);
age_ERR = 2*mean(std(radar.age, [], 3), 2);

figure
hold on
h1 = plot(core.age, core.depth, 'b', 'LineWidth', 2);
h2 = plot(age_mean, radar.depth_interp, 'r', 'LineWidth', 2);
h3 = plot(age_mean+age_ERR, radar.depth_interp, 'r--', 'LineWidth', 0.5);
h4 = plot(age_mean-age_ERR, radar.depth_interp, 'r--', 'LineWidth', 0.5);
xlabel('Calendar Year')
ylabel('Depth (m)')
xlim([min([min(core.age) min(age_mean-age_ERR)]) max([max(core.age) max(age_mean)])])
legend([h1 h2], 'Core age (manual)', 'Radar age (automated)', 'Location', 'nw')
set(gca, 'Ydir', 'reverse', 'FontSize', 18)
hold off
