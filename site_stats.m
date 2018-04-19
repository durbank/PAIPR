% This is a script for investigating statistical relationships for radar
% from individual firn core sites

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/Data/';
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

% Define number of MC realizations
Ndraw = 100;

% Import firn core data
[cores] = import_cores(strcat(data_path, 'Ice-cores/SEAT_cores/', ...
    'DGK_core_data.xlsx'), Ndraw);

%%

inputs = {'SEAT10_1' 'SEAT10_3'};

for i = 1:length(inputs)
    file = strcat(data_path, 'radar/SEAT_Traverses/core-site_tests/', ...
        inputs{i}, '_radar.mat');
    radar.(inputs{i}) = load(file);
end

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
