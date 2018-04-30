


% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'F:/Research/Antarctica/Data/';
        addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Addons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Number of Monte Carlo simulations
Ndraw = 100;

% Import firn core data
[cores] = import_cores(strcat(data_path, 'Ice-cores/SEAT_cores/', ...
    'DGK_core_data.xlsx'), Ndraw);

%% Define radar files to import/process

radar_dir = strcat(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband/', ...
    'ProcessedSEAT2010/transectSEAT10_1_2/');

% List all files matching 'wild' within radar directory
wild = 'layers*';
files = dir(strcat(radar_dir, wild));

i = randi(length(files));
file = strcat(radar_dir, files(i).name);

% Path to full SEAT transect
file = strcat(data_path, 'radar/SEAT_Traverses/core-site_tests/', ...
    'layers_ku_band_SEAT10_5.mat');

% % Path of the OIB file to process
% % SEAT10_4
% file = strcat(data_path, 'IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_272.nc');
% file = strcat(data_path, 'IceBridge/Snow Radar/2016/IRSNO1B_20161109_02_381.nc');
% file = strcat(data_path, 'IceBridge/Kuband/2016/IRKUB1B_20161109_02_381.nc');
% % SEAT10_5
% file = strcat(data_path, 'IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_257.nc');
% % SEAT10_6
% file = strcat(data_path, 'IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_242.nc');

%%

% for i = 1:length(files)

% Calculate radar ages and associated other data
[radar] = radar_age(file, cores, Ndraw);

% Calculate annual accumulation rates from data
[radar0] = calc_SWE(radar, Ndraw);

radar = radar0;

% Remove first/last 10 traces in radar (addresses some edge effect problems
% present in many data sets
edge = 25;

radar.Easting = radar0.Easting(edge:end-edge);
radar.Northing = radar0.Northing(edge:end-edge);
radar.dist = radar0.dist(edge:end-edge);
radar.data_stack = radar0.data_stack(:,edge:end-edge);
radar.rho_coeff = radar0.rho_coeff(:,edge:end-edge);
radar.rho_var = radar0.rho_var(:,edge:end-edge);
radar.data_smooth = radar0.data_smooth(:,edge:end-edge);
radar.layer_vals = radar0.layer_vals(:,edge:end-edge);
radar.likelihood = radar0.likelihood(:,edge:end-edge);
radar.ages = radar0.ages(:,edge:end-edge,:);
radar.SMB_yr = radar0.SMB_yr(edge:end-edge);
radar.SMB = radar0.SMB(edge:end-edge);

%%

% Calculate mean accumulation rate and std at each location
SMB_mean = cellfun(@mean, radar.SMB, 'UniformOutput', 0);
SMB_std = cellfun(@std, radar.SMB, 'UniformOutput', 0);

trend = cell(1, length(radar.SMB));
p_val = cell(1, length(radar.SMB));
for i = 1:length(radar.SMB)
    
    trend_i = zeros(1, Ndraw);
    p_val_i = zeros(1, Ndraw);
    for j = 1:Ndraw
        % Regression using regress function (average of MC simulations)
        [b, ~, ~, ~, stats] = regress(radar.SMB{i}(:,j), ...
            [ones(length(radar.SMB_yr{i}), 1) radar.SMB_yr{i}]);
        trend_i(j) = b(2);
        p_val_i(j) = stats(3);
    end
    trend{i} = trend_i;
    p_val{i} = p_val_i;
end

trend_mean = cellfun(@mean, trend);
trend_std = cellfun(@std, trend);

% Find indices of significant trends (at 95% confidence level)
p_star = cellfun(@(x) find(x<=0.05), p_val, 'UniformOutput', 0);

% Calculate percentage of MC realizations with significant trends
p_ratio = cellfun(@(sig, all) length(sig)/length(all), p_star, p_val);

% Regression of SMB trend on mean SMB (all traces)
[trend_b, trend_int, ~, ~, trend_stats] = regress(trend_mean', ...
    [ones(length(trend_mean), 1) cellfun(@mean, SMB_mean)']);
[p_b, S_b] = polyfit(cellfun(@mean, SMB_mean), trend_mean, 1);
[trend_Y, trend_D] = polyconf(p_b, cellfun(@mean, SMB_mean), S_b);

%% Diagnostic plots for single random trace

% Trace idx to investigate
i = randi(length(radar.SMB));

% Find the nearest cores to the radar data (for comparison plots)
[~, cores_near_idx] = sort(pdist2([radar.Easting(i) radar.Northing(i)], ...
    [cores.Easting' cores.Northing'], 'Euclidean'));
core_near1 = cores.(cores.name{cores_near_idx(1)});
core_near2 = cores.(cores.name{cores_near_idx(2)});
core_near3 = cores.(cores.name{cores_near_idx(3)});

% Calculate the mean age-depth scale and std for radar trace i
age_mean = mean(squeeze(radar.ages(:,i,:)), 2);
age_std = std(squeeze(radar.ages(:,i,:)), [], 2);

% Plot full radargram
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

% Age-depth scale comparison between radar trace and nearest cores
figure
hold on
h1 = plot(core_near1.depth, mean(core_near1.ages, 2), 'b', 'LineWidth', 2);
plot(core_near1.depth, mean(core_near1.ages, 2) + 2*std(core_near1.ages, [], 2), 'b--')
plot(core_near1.depth, mean(core_near1.ages, 2) - 2*std(core_near1.ages, [], 2), 'b--')
h2 = plot(core_near2.depth, core_near2.age, 'c', 'LineWidth', 2);
h3 = plot(core_near3.depth, core_near3.age, 'c--', 'LineWidth', 1);
h4 = plot(radar.depth, age_mean, 'r', 'LineWidth', 2);
plot(radar.depth, age_mean + 2*age_std, 'r--', 'LineWidth', 0.5)
plot(radar.depth, age_mean - 2*age_std, 'r--', 'LineWidth', 0.5)
ylabel('Calendar Year')
xlabel('Depth (m)')
legend([h1 h2 h3 h4], 'Nearest core age (manual)', '2nd nearest core', ...
    '3rd nearest core', 'Radar age (automated)', 'Location', 'ne')
set(gca, 'FontSize', 10)
hold off

% Compare annual accumulation between ith radar trace and nearest 2 cores
% (along with estimated uncertainties for each)
figure
hold on
h1 = plot(core_near1.SMB_yr, mean(core_near1.SMB, 2), 'b', 'LineWidth', 2);
plot(core_near1.SMB_yr, mean(core_near1.SMB, 2) + 2*std(core_near1.SMB, [], 2), 'b--')
plot(core_near1.SMB_yr, mean(core_near1.SMB, 2) - 2*std(core_near1.SMB, [], 2), 'b--')

h2 = plot(core_near2.SMB_yr, mean(core_near2.SMB, 2), 'c', 'LineWidth', 2);
plot(core_near2.SMB_yr, mean(core_near2.SMB, 2) + 2*std(core_near2.SMB, [], 2), 'c--')
plot(core_near2.SMB_yr, mean(core_near2.SMB, 2) - 2*std(core_near2.SMB, [], 2), 'c--')

h3 = plot(radar.SMB_yr{i}, mean(radar.SMB{i}, 2), 'r', 'LineWidth', 2);
plot(radar.SMB_yr{i}, mean(radar.SMB{i}, 2) + 2*std(radar.SMB{i}, [], 2), 'r--')
plot(radar.SMB_yr{i}, mean(radar.SMB{i}, 2) - 2*std(radar.SMB{i}, [], 2), 'r--')
legend([h1 h2 h3], 'Nearest firn core', '2nd nearest core', 'Ku radar')
xlabel('Calendar Year')
ylabel('Annual accumulation (mm w.e.)')
hold off

% Compare linear trend in SMB for ith trace and nearest 3 cores (with
% uncertainty)



%% Diagnostic plots for bulk radar file

% Mean SMB across entire radargram (mean of all realizations for each trace)
figure
scatter(radar.Easting, radar.Northing, 30, cellfun(@mean, SMB_mean), 'filled')
hcb = colorbar;
ylabel(hcb, 'Mean annual SMB (mm/a')
colormap(cool)

% Mean trend in SMB across entire radargram
figure
scatter(radar.Easting, radar.Northing, 30, trend_mean, 'filled')
hcb = colorbar;
ylabel(hcb, 'Trend in SMB (mm/a)')
colormap(cool)

% Mean SMB vs mean trend
figure
hold on
plot(cellfun(@mean, SMB_mean), trend_mean, 'b.')
h1 = plot(cellfun(@mean, SMB_mean), trend_Y, 'b', 'LineWidth', 2);
plot(cellfun(@mean, SMB_mean), trend_Y + trend_D, 'b--')
plot(cellfun(@mean, SMB_mean), trend_Y - trend_D, 'b--')
xlabel('Mean annual accumulation (mm w.e.)')
ylabel('Annual accumulation trend (mm/a)')
hold off