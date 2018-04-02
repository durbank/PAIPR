


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

% Find the nearest core to the radar data (for comparison plots)
trace_idx = round(size(radar.data_smooth, 2)/2);
[~, core_near_idx] = min(pdist2([radar.Easting(trace_idx) ...
    radar.Northing(trace_idx)], [cores.Easting' cores.Northing'], ...
    'Euclidean'));
core_near = cores.(cores.name{core_near_idx});

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


% Calculate SWE accumulation at each depth interval in the weighted
% composite core
core_accum_dt = 0.02*(1000*core_near.rho);

% Find indices of integer ages within core age profile
yr_top = floor(core_near.age(2));
yr_end = ceil(core_near.age(end));
core_yr = (yr_top:-1:yr_end)';
core_yr_idx = logical([1; diff(floor(core_near.age))]);
yr_loc = find(core_yr_idx);

core_accum = zeros(length(core_yr), Ndraw);
for j = 1:Ndraw
    
    % Add noise to integer age locations due to uncertainty in exact point 
    % in time of the accumulation peak, using a std dev of 1 month
    yr_loc_j = yr_loc;
    yr_loc_j(2:end-1) = yr_loc(2:end-1) + ...
        round(1*(mean(diff(yr_loc))/12)*randn(length(yr_loc)-2, 1));
    loc_idx = yr_loc_j<1;
    yr_loc_j(loc_idx) = yr_loc(loc_idx);
    
    % Integrate accumulation at each depth point for each whole year in
    % firn core
    core_accum_j = zeros(length(core_yr), 1);
    for n = 1:length(core_yr)
        core_accum_j(n) = sum(core_accum_dt(yr_loc_j(n)+1:yr_loc_j(n+1)));
    end
    
    % Output accumulatio results to preallocated array
    core_accum(:,j) = core_accum_j;
end
core_near.SMB_yr = core_yr;
core_near.SMB = core_accum;


accum_mean = mean(cell2mat(radar.SMB(:,i)), 2);
accum_ERR = 2*std(cell2mat(radar.SMB(:,i)), [], 2);

figure
hold on
h1 = plot(core_near.SMB_yr, mean(core_near.SMB, 2), 'b', 'LineWidth', 2);
plot(core_near.SMB_yr, mean(core_near.SMB, 2) + 2*std(core_near.SMB, [], 2), 'b--')
plot(core_near.SMB_yr, mean(core_near.SMB, 2) - 2*std(core_near.SMB, [], 2), 'b--')
h2 = plot(radar.SMB_yr{i}, accum_mean, 'r', 'LineWidth', 2);
plot(radar.SMB_yr{i}, accum_mean + accum_ERR, 'r--')
plot(radar.SMB_yr{i}, accum_mean - accum_ERR, 'r--')
legend([h1 h2], 'Firn core', 'Ku radar')
xlabel('Calendar Year')
ylabel('Annual accumulation (mm w.e.)')
hold off
