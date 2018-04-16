


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

% % Path to full SEAT transect
% file = strcat(data_path, 'OUTPUT\SEAT2010_transects\', ...
%     'layers_ku_band_transectSEAT10_5_6.mat');
% 
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
radar.age = radar0.age(:,edge:end-edge,:);
radar.SMB_yr = radar0.SMB_yr(edge:end-edge);
radar.SMB = radar0.SMB(edge:end-edge);

%%

% Calculate mean accumulation rate and std at each location
SMB_mean = cellfun(@mean, radar.SMB, 'UniformOutput', 0);
SMB_std = cellfun(@std, radar.SMB, 'UniformOutput', 0);

figure
scatter(radar.Easting, radar.Northing, 30, cellfun(@mean, SMB_mean), 'filled')
hcb = colorbar;
ylabel(hcb, 'Mean annual SMB (mm/a')
colormap(cool)

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

% Plot trends
figure
scatter(radar.Easting, radar.Northing, 30, trend_mean, 'filled')
hcb = colorbar;
ylabel(hcb, 'Trend in SMB (mm/a)')
colormap(cool)

%%

% Regression of SMB trend on mean SMB (all traces)
[trend_b, trend_int, ~, ~, trend_stats] = regress(trend_mean', ...
    [ones(length(trend_mean), 1) cellfun(@mean, SMB_mean)']);
[p_b, S_b] = polyfit(cellfun(@mean, SMB_mean), trend_mean, 1);
[trend_Y, trend_D] = polyconf(p_b, cellfun(@mean, SMB_mean), S_b);

% % Regression of SMB trend of mean SMB (traces with significant trends)
% [star_trend, star_int, ~, ~, star_stats] = regress(trend(p_star)', ...
%     [ones(length(trend(p_star)), 1) SMB_mean(p_star)']);
% [coeff_star, S_star] = polyfit(SMB_mean(p_star), trend(p_star), 1);
% [star_Y, star_D] = polyconf(coeff_star, SMB_mean(p_star), S_star);

% Plot mean SMB vs trend (for all trends, and significant trends
figure
hold on
plot(cellfun(@mean, SMB_mean), trend_mean, 'b.')
h1 = plot(cellfun(@mean, SMB_mean), trend_Y, 'b', 'LineWidth', 2);
plot(cellfun(@mean, SMB_mean), trend_Y + trend_D, 'b--')
plot(cellfun(@mean, SMB_mean), trend_Y - trend_D, 'b--')
% plot(SMB_mean(p_star), trend(p_star), 'r.')
% h2 = plot(SMB_mean(p_star), star_Y, 'r', 'LineWidth', 2);
% plot(SMB_mean(p_star), star_Y + star_D, 'r--')
% plot(SMB_mean(p_star), star_Y - star_D, 'r--')
% legend([h1 h2], 'All traces', 'Significant trends')
xlabel('Mean annual accumulation (mm w.e.)')
ylabel('Annual accumulation trend (mm/a)')
hold off

%%

%  figure
% yyaxis left
% plot(SMB_mean, 100*trend./SMB_mean, 'bo')
% hold on
% yyaxis right
% plot(SMB_mean, trend, 'ro')
% hold off



% [p, S] = cellfun(@(year, SMB) polyfit(year, mean(SMB, 2), 1), radar.SMB_yr, radar.SMB, 'UniformOutput', 0);
% Y = cellfun(@(fit, X) polyconf(fit, X), p, radar.SMB_yr, 'UniformOutput', 0);
% idx = 1:10:length(radar.SMB);
% 
% figure
% hold on
% for i = idx
%     plot(radar.SMB_yr{i}, Y{i})
% end
% hold off

% figure
% hold on
% plot(cores.SEAT10_4.SMB_yr, mean(cores.SEAT10_4.SMB, 2), 'k')
% 
% for i = idx
%     plot(radar.SMB_yr{i}, mean(radar.SMB{i}, 2))
% end
% hold off
