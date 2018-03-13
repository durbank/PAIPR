% Directory to radar files of interest

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

% Import firn core data
[cores] = import_cores(strcat(data_path, ['SEAT_cores' filesep ...
    'DGK_core_data.xlsx']));

%% Define radar file to import/process

% radar_dir = strcat(data_path, ['SEAT_Traverses' filesep 'SEAT2010Kuband'...
%     filesep 'ProcessedSEAT2010' filesep 'grid_SEAT10_4' filesep]);
% 
% % List all files matching 'wild' within radar directory
% wild = 'layers*';
% files = dir(strcat(radar_dir, wild));
% 
% % Select individual radar file for analysis
% i = randi(length(files));
% file = strcat([files(i).folder filesep], files(i).name);

% Path of the OIB file to process
file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Kuband/2016/IRKUB1B_20161109_02_381.nc';

%%
% Number of simulations to perform on age-depth Monte Carlo
Ndraw = 100;

% Calculate radar ages and associated other data
[radar, core] = age_MC(file, cores, Ndraw);

% Calculate annual accumulation rates from data
[SMB] = calc_SWE(radar, core);

% Diagnostic figure

trace_near = round(size(SMB.radar_accum, 2)/2);
figure
hold on
h1 = plot(SMB.radar_yr, SMB.radar_accum(:,trace_near), 'm--', 'LineWidth', 0.25);
h2 = plot(SMB.core_yr, SMB.core_accum, 'k--');
h3 = plot(SMB.core_yr, movmean(SMB.core_accum, 3), 'b', 'LineWidth', 2);
h4 = plot(SMB.radar_yr, mean(SMB.radar_accum, 2), 'r', 'LineWidth', 2);
plot(SMB.radar_yr, mean(SMB.radar_accum, 2) + std(SMB.radar_accum, [], 2), 'r--', 'LineWidth', 2)
plot(SMB.radar_yr, mean(SMB.radar_accum, 2) - std(SMB.radar_accum, [], 2), 'r--', 'LineWidth', 2)
legend([h1 h2 h3 h4], 'Nearest trace', 'Core accum', ...
    'Core 3-yr moving mean', 'Mean radar accum')
xlabel('Year C.E.')
ylabel('Accumulation [mm w.e./a]')
xlim([min([min(SMB.radar_yr) min(SMB.core_yr)]) max([max(SMB.radar_yr) max(SMB.core_yr)])])
hold off
