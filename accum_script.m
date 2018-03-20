


% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/WAIS Variability/';
        addon_path = 'E:/Research/Antarctica/WAIS Variability/Addons/';
    case false
        data_path = '/Volumes/WARP/Research/Antarctica/WAIS Variability/';
        addon_path = '/Users/Durbank/Documents/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Import firn core data
[cores] = import_cores(strcat(data_path, ['SEAT_cores' filesep ...
    'DGK_core_data.xlsx']));

%% Define radar file to import/process

radar_dir = strcat(data_path, ['SEAT_Traverses' filesep 'SEAT2010Kuband'...
    filesep 'ProcessedSEAT2010' filesep 'grid_SEAT10_4' filesep]);

% List all files matching 'wild' within radar directory
wild = 'layers*';
files = dir(strcat(radar_dir, wild));

% Select individual radar file for analysis
i = randi(length(files));
file = strcat([files(i).folder filesep], files(i).name);

% Path of the OIB file to process
% SEAT10_4
% file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_272.nc';
% file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2016/IRSNO1B_20161109_02_381.nc';
% file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Kuband/2016/IRKUB1B_20161109_02_381.nc';
% SEAT10_5
% file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_257.nc';
% SEAT10_6
% file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_242.nc';

%%
% Number of simulations to perform on age-depth Monte Carlo
Ndraw = 100;

% Calculate radar ages and associated other data
[radar, core] = radar_age(file, cores, Ndraw);

% Calculate annual accumulation rates from data
[radar, core] = calc_SWE(radar, core, Ndraw);

%% Diagnostic figure

% Select random radar trace for comparison plots
i = randi(size(radar.data_smooth, 2));

% Plot radargram
figure('Position', [200 200 1500 800])
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
colorbar
xlabel('Distance along profile (m)')
ylabel('Depth (m)')
hold on
plot([radar.dist(i) radar.dist(i)], [0 radar.depth(end)], 'r', 'LineWidth', 2)
xlim([0 radar.dist(end)])
ylim([0 radar.depth(end)])
set(gca, 'Ydir', 'reverse', 'FontSize', 18)
hold off

% Compare radar signal at random trace to annual horizion selection at that
% trace
% layer_idx = logical([diff(floor(squeeze(radar.age(:,i,100)))); 0]);
layer_idx = logical([diff(floor(squeeze(radar.age(:,i,randi(Ndraw))))); 0]);
figure('Position', [50 50 500 1200])
hold on
plot(radar.data_smooth(:,i), radar.depth, 'r', 'LineWidth', 1.5)
scatter(radar.data_smooth(layer_idx,i), radar.depth(layer_idx), 'b<', 'filled')
xlim([-1.25 2])
xlabel('Radar Z-score')
ylabel('Depth (m)')
set(gca, 'Ydir', 'reverse')
hold off


figure
plot(radar.depth, radar.layer_vals(:,i))


age_mean = mean(squeeze(radar.age(:,i,:)), 2);
age_ERR = 2*std(squeeze(radar.age(:,i,:)), [], 2);

figure
hold on
h1 = plot(core.depth, core.age, 'b', 'LineWidth', 2);
h2 = plot(radar.depth, age_mean, 'r', 'LineWidth', 2);
plot(radar.depth, age_mean + age_ERR, 'r--', 'LineWidth', 0.5)
plot(radar.depth, age_mean - age_ERR, 'r--', 'LineWidth', 0.5)
ylabel('Calendar Year')
xlabel('Depth (m)')
ylim([min([min(core.age) min(age_mean-age_ERR)]) max([max(core.age) max(age_mean)])])
legend([h1 h2], 'Core age (manual)', 'Radar age (automated)', 'Location', 'ne')
set(gca, 'FontSize', 10)
hold off

figure
hold on
plot(core.SMB_yr, core.SMB, 'b')
plot(radar.SMB_yr, radar.SMB(:,i), 'r')
hold off


% trace_near = round(size(SMB.radar_accum, 2)/2);
% figure
% hold on
% h1 = plot(SMB.radar_yr, SMB.radar_accum(:,trace_near), 'm--', 'LineWidth', 0.25);
% h2 = plot(SMB.core_yr, SMB.core_accum, 'k--');
% h3 = plot(SMB.core_yr, movmean(SMB.core_accum, 3), 'b', 'LineWidth', 2);
% h4 = plot(SMB.radar_yr, mean(SMB.radar_accum, 2), 'r', 'LineWidth', 2);
% plot(SMB.radar_yr, mean(SMB.radar_accum, 2) + std(SMB.radar_accum, [], 2), 'r--', 'LineWidth', 2)
% plot(SMB.radar_yr, mean(SMB.radar_accum, 2) - std(SMB.radar_accum, [], 2), 'r--', 'LineWidth', 2)
% legend([h1 h2 h3 h4], 'Nearest trace', 'Core accum', ...
%     'Core 3-yr moving mean', 'Mean radar accum')
% xlabel('Year C.E.')
% ylabel('Accumulation [mm w.e./a]')
% xlim([min([min(SMB.radar_yr) min(SMB.core_yr)]) max([max(SMB.radar_yr) max(SMB.core_yr)])])
% hold off
