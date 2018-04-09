


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

%% Define radar files to import/process

radar_dir = strcat(data_path, ['SEAT_Traverses' filesep 'SEAT2010Kuband'...
    filesep 'ProcessedSEAT2010' filesep 'transectSEAT10_1_2' filesep]);

% List all files matching 'wild' within radar directory
wild = 'layers*';
files = dir(strcat(radar_dir, wild));

%%
% Number of simulations to perform on age-depth Monte Carlo
Ndraw = 100;

% for i = 1:length(files)

i = 2;
    file = strcat(radar_dir, files(i).name);

% Calculate radar ages and associated other data
[radar, core] = radar_age(file, cores, Ndraw);

% Calculate annual accumulation rates from data
[radar, core] = calc_SWE(radar, core, Ndraw);

% Calculate mean accumulation rate and std at each location
SMB_mean = cellfun(@(x) mean(mean(x)), radar.SMB);
SMB_std = cellfun(@(x) mean(std(x)), radar.SMB);
% SMB_std = cellfun(@(x) mean(std(x)/sqrt(length(x))), radar.SMB);
grid_mean = mean(SMB_mean);
grid_ERR = 1.96*std(SMB_mean)/sqrt(length(SMB_mean));

% % Calculate linear trend in accumulation rate and uncertainty at each
% % location
% [P, err] = cellfun(@(x,y) polyfit(x,mean(y, 2),1), radar.SMB_yr, ...
%     radar.SMB, 'UniformOutput', 0);
% [trendline, trend_std] = cellfun(@(x,y,z) polyval(x, y, z), ...
%     P, radar.SMB_yr, err, 'UniformOutput', 0);

% Regression using regress function
[b, bint, r] = cellfun(@(SMB, year) regress(mean(SMB, 2), year), radar.SMB, ...
    radar.SMB_yr, 'UniformOutput', 0);

