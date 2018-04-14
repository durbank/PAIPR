


% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/WAIS Variability/';
        addon_path = 'E:/Research/Antarctica/WAIS Variability/Addons/';
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
[cores] = import_cores(strcat(data_path, ['Ice-cores' filesep 'SEAT_cores' ...
    filesep 'DGK_core_data.xlsx']), Ndraw);

%% Define radar files to import/process

radar_dir = strcat(data_path, ['radar' filesep 'SEAT_Traverses' filesep ...
    'SEAT2010Kuband' filesep 'ProcessedSEAT2010' filesep ...
    'transectSEAT10_5_6' filesep]);

% List all files matching 'wild' within radar directory
wild = 'layers*';
files = dir(strcat(radar_dir, wild));

i = randi(length(files));
file = strcat(radar_dir, files(i).name);

% Path to full SEAT transect
% file = strcat(data_path, 'OUTPUT\SEAT2010_transects\layers_ku_band_transectSEAT10_5_6.mat');
% file = 'E:\Research\Antarctica\Data\OUTPUT\SEAT2010_transects\layers_ku_band_transectSEAT10_5_6.mat';

% Path of the OIB file to process
% SEAT10_4
file = strcat(data_path, 'IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_272.nc');
% file = 'E:\Research/Antarctica/Data/IceBridge/Snow Radar/2016/IRSNO1B_20161109_02_381.nc';
% file = 'E:\Research/Antarctica/Data/IceBridge/Kuband/2016/IRKUB1B_20161109_02_381.nc';
% SEAT10_5
% file = 'E:\Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_257.nc';
% SEAT10_6
% file = 'E:\Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_242.nc';

%%

% for i = 1:length(files)

% Calculate radar ages and associated other data
[radar] = radar_age(file, cores, Ndraw);

% Calculate annual accumulation rates from data
[radar0] = calc_SWE(radar, Ndraw);

radar = radar0;

% Remove first/last 10 traces in radar (addresses some edge effect problems
% present in many data sets
radar.Easting = radar0.Easting(10:end-10);
radar.Northing = radar0.Northing(10:end-10);
radar.dist = radar0.dist(10:end-10);
radar.data_stack = radar0.data_stack(:,10:end-10);
radar.rho_coeff = radar0.rho_coeff(:,10:end-10);
radar.rho_var = radar0.rho_var(:,10:end-10);
radar.data_smooth = radar0.data_smooth(:,10:end-10);
radar.layer_vals = radar0.layer_vals(:,10:end-10);
radar.likelihood = radar0.likelihood(:,10:end-10);
radar.age = radar0.age(:,10:end-10,:);
radar.SMB_yr = radar0.SMB_yr(10:end-10);
radar.SMB = radar0.SMB(10:end-10);

%%

% Calculate mean accumulation rate and std at each location
SMB_mean = cellfun(@(x) mean(mean(x)), radar.SMB);
SMB_std = cellfun(@(x) mean(std(x)), radar.SMB);

figure
scatter(radar.Easting, radar.Northing, 30, SMB_mean, 'filled')
hcb = colorbar;
ylabel(hcb, 'Mean annual SMB (mm/a')
colormap(cool)

% Regression using regress function (average of MC simulations)
[b, ~, r, ~, stats] = cellfun(@(SMB, year) regress(mean(SMB, 2), ...
    [ones(length(year),1) year]), radar.SMB, radar.SMB_yr,...
    'UniformOutput', 0);
% ***REMOVES FIRST YEAR DUE TO ISSUES WITH HIGHLY REPRESSED ACCUMULATION 
% RATES FOR TOP SURFACE YEAR ESTIMATES***
% [b, ~, r, ~, stats] = cellfun(@(SMB, year) regress(mean(SMB, 2), ...
%     [ones(length(year),1) year]), radar.SMB, radar.SMB_yr,...
%     'UniformOutput', 0);

% Extract trends from b variable
trend = cellfun(@(x) x(2), b);

% Extract p-values from stats
p_val = cellfun(@(x) x(3), stats);

% Find indices of significant trends (at 95% confidence level)
p_star = find(p_val<=0.05);

% Calculate percentage of traces with significant trend
p_ratio = length(p_star)/length(p_val);

% % Create colormap that is red for negative, blue for positive,
% % and white in the middle
% redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
% blueColorMap = [zeros(1, 132), linspace(0, 1, 124)];
% colorMap = [redColorMap; zeros(1, 256); blueColorMap]';
% colorMap( ~any(colorMap,2), : ) = 1; 

% Plot trends
figure
scatter(radar.Easting, radar.Northing, 30, trend, 'filled')
hcb = colorbar;
ylabel(hcb, 'Trend in SMB (mm/a')
colormap(cool)

% Regression of SMB trend on mean SMB (all traces)
[trend_b, trend_int, ~, ~, trend_stats] = regress(trend', ...
    [ones(length(trend), 1) SMB_mean']);
[p_b, S_b] = polyfit(SMB_mean, trend, 1);
[trend_Y, trend_D] = polyconf(p_b, SMB_mean, S_b);

% Regression of SMB trend of mean SMB (traces with significant trends)
[star_trend, star_int, ~, ~, star_stats] = regress(trend(p_star)', ...
    [ones(length(trend(p_star)), 1) SMB_mean(p_star)']);
[coeff_star, S_star] = polyfit(SMB_mean(p_star), trend(p_star), 1);
[star_Y, star_D] = polyconf(coeff_star, SMB_mean(p_star), S_star);

% Plot mean SMB vs trend (for all trends, and significant trends
figure
hold on
plot(SMB_mean, trend, 'b.')
h1 = plot(SMB_mean, trend_Y, 'b', 'LineWidth', 2);
plot(SMB_mean, trend_Y + trend_D, 'b--')
plot(SMB_mean, trend_Y - trend_D, 'b--')
plot(SMB_mean(p_star), trend(p_star), 'r.')
h2 = plot(SMB_mean(p_star), star_Y, 'r', 'LineWidth', 2);
plot(SMB_mean(p_star), star_Y + star_D, 'r--')
plot(SMB_mean(p_star), star_Y - star_D, 'r--')
legend([h1 h2], 'All traces', 'Significant trends')
xlabel('Mean annual accumulation (mm w.e.)')
ylabel('Annual accumulation trend (mm/a)')
hold off

 figure
yyaxis left
plot(SMB_mean, 100*trend./SMB_mean, 'bo')
hold on
yyaxis right
plot(SMB_mean, trend, 'ro')
hold off



[p, S] = cellfun(@(year, SMB) polyfit(year, mean(SMB, 2), 1), radar.SMB_yr, radar.SMB, 'UniformOutput', 0);
Y = cellfun(@(fit, X) polyconf(fit, X), p, radar.SMB_yr, 'UniformOutput', 0);
idx = 1:10:length(radar.SMB);

figure
hold on
for i = idx
    plot(radar.SMB_yr{i}, Y{i})
end