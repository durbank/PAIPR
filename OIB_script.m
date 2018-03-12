% Script for testing the import and processing of OIB data with my codes

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

addpath cresis-L1B-matlab-readers/

% Files IRSNO1B_20111109_02_211 through ...02_227 roughly follows the
% SEAT10-1 to SEAT10-6 transect. The same holds true for the following:
%   IRSNO1B_20161109_02_320 through IRSNO1B_20161109_02_336
%   Ku-band: IRKUB1B_20161109_02_320 through IRKUB1B_20161109_02_336

% Files IRSNO1B_20111109_02_242 through ...02_272 exactly follows the 
% SEAT10-4 - SEAT10-6 transect. The same holds true for the following:
%   IRSNO1B_20161109_02_350 through IRSNO1B_20161109_02_381
%   Ku-band: IRKUB1B_20161109_02_350 through IRKUB1B_20161109_02_381

% Path of the OIB file to process
file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow radar/2011/IRSNO1B_20111109_02_211.nc';

% Number of simulations to perform on age-depth Monte Carlo
Ndraw = 100;

%%

% Import firn core data
[cores] = import_cores(strcat(data_path, ['SEAT_cores' filesep ...
    'DGK_core_data.xlsx']));

% Calculate radar ages and associated other data
[radar, core] = age_MC(file, cores, Ndraw);

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

% Plot core vs mean radar ages (and uncertainty) for random trace
age_mean = mean(squeeze(radar.age(:,i,:)), 2);
age_ERR = 2*std(squeeze(radar.age(:,i,:)), [], 2);

figure
hold on
h1 = plot(core.depth, core.age, 'b', 'LineWidth', 2);
h2 = plot(radar.depth_interp, age_mean, 'r', 'LineWidth', 2);
h3 = plot(radar.depth_interp, age_mean+age_ERR, 'r--', 'LineWidth', 0.5);
h4 = plot(radar.depth_interp, age_mean-age_ERR, 'r--', 'LineWidth', 0.5);
xlabel('Depth (m)')
ylabel('Calendar Year')
ylim([min([min(core.age) min(age_mean-age_ERR)]) max([max(core.age) max(age_mean)])])
legend([h1 h2], 'Core age (manual)', 'Radar age (automated)', 'Location', 'ne')
set(gca, 'FontSize', 18)
hold off


