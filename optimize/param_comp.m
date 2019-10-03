% Script to compare the results of logistic regression parameter estimation
% between the different validation sites

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/Data/';
        addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile(addon_path, 'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add PAIPR-core functions to path
PAIPR_path = fullfile(cd,'..', 'PAIPR-core');
addpath(genpath(PAIPR_path))

%%

input_dir = fullfile(data_path, "PAIPR-results/vWIP",...
    "coeffs_branch/manual_layers");

% SEAT10_4 = load(fullfile(input_dir, "SEAT10_4", ...
%     "params_scaled.mat"));
% SEAT10_5 = load(fullfile(input_dir, "SEAT10_5", ...
%     "params_scaled.mat"));
% SEAT10_6 = load(fullfile(input_dir, "SEAT10_6", ...
%     "params_scaled.mat"));
SEAT10_4 = load(fullfile(input_dir, "SEAT10_4", ...
    "params_output.mat"));
SEAT10_5 = load(fullfile(input_dir, "SEAT10_5", ...
    "params_output.mat"));
SEAT10_6 = load(fullfile(input_dir, "SEAT10_6", ...
    "params_output.mat"));

%%
figure
title("Rate parameter (r)")
hold on
ksdensity(SEAT10_4.r_params)
ksdensity(SEAT10_5.r_params)
ksdensity(SEAT10_6.r_params)
legend("2010-4", "2010-5", "2010-6")
% xlim([-6.5 -3.5])
hold off

% ksdensity()
% findpeaks()

sprintf("Median r parameters: %0.2e | %0.2e | %0.2e", ...
    median(SEAT10_4.r_params), median(SEAT10_5.r_params), ...
    median(SEAT10_6.r_params))

% figure
% title("Midpoint parameter (k)")
% hold on
% ksdensity(SEAT10_4.k_params)
% ksdensity(SEAT10_5.k_params)
% ksdensity(SEAT10_6.k_params)
% legend("2010-4", "2010-5", "2010-6")
% % xlim([-2 3])
% hold off
% sprintf("Median k parameters: %0.2f | %0.2f | %0.2f", ...
%     median(SEAT10_4.k_params), median(SEAT10_5.k_params), ...
%     median(SEAT10_6.k_params))

figure
title("SSE")
hold on
ksdensity(SEAT10_4.SSE)
ksdensity(SEAT10_5.SSE)
ksdensity(SEAT10_6.SSE)
legend("2010-4", "2010-5", "2010-6")
% xlim([0 50])
hold off
sprintf("Median SSE values: %0.2f | %0.2f | %0.2f", ...
    median(SEAT10_4.SSE), median(SEAT10_5.SSE), ...
    median(SEAT10_6.SSE))

%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

r = mean([median(SEAT10_4.r_params), median(SEAT10_5.r_params), ...
    median(SEAT10_6.r_params)]);
k = mean([median(SEAT10_4.k_params), median(SEAT10_5.k_params), ...
    median(SEAT10_6.k_params)]);

% % max bounds
% r = -2.33e-4;
% k = 3.75;
% % min bounds
% r = -3.67e-4;
% k = 2.25;

% Calculate ages for validation sites
radar4 = load(fullfile(input_dir, "SEAT10_4", "PAIPR_out.mat"));
[radar4] = radar_age(radar4, r, k, 100);
radar5 = load(fullfile(input_dir, "SEAT10_5", "PAIPR_out.mat"));
[radar5] = radar_age(radar5, r, k, 100);
radar6 = load(fullfile(input_dir, "SEAT10_6", "PAIPR_out.mat"));
[radar6] = radar_age(radar6, r, k, 100);


X = 0:(25*length(radar4.Easting)*max(radar4.peaks(:)));
% X = 0:0.01:3.5;
Y = 1./(1+exp(r*X + k));
Y_ub = 1./(1+exp(-2.33e-4*X + 3.75));
Y_lb = 1./(1+exp(-3.67e-4*X + 2.25));

figure
hold on
plot(radar4.DB(radar4.DB>0), radar4.likelihood(radar4.DB>0), 'o')
plot(radar5.DB(radar5.DB>0), radar5.likelihood(radar5.DB>0), 'o')
plot(radar6.DB(radar6.DB>0), radar6.likelihood(radar6.DB>0), 'o')
plot(X,Y, 'b--')
plot(X,Y_lb, 'r--')
plot(X,Y_ub, 'r--')
plot(radar4.DB(radar4.DB>0), radar4.likelihood(radar4.DB>0), 'o')


%%

% Age plots for SEAT2010-4
trace = round(length(radar4.Easting)/2);
age4 = mean(squeeze(radar4.ages(:,trace,:)),2);
std4 = std(squeeze(radar4.ages(:,trace,:)),[],2);
figure
hold on
plot(radar4.depth, age4, 'r', 'LineWidth', 2)
plot(radar4.depth, age4+std4, 'r--')
plot(radar4.depth, age4-std4, 'r--')
plot(cores.SEAT10_4.depth, mean(cores.SEAT10_4.ages,2),...
    'b', 'LineWidth', 2)
plot(cores.SEAT10_4.depth, mean(cores.SEAT10_4.ages,2) + ...
    std(cores.SEAT10_4.ages,[],2), 'b--')
plot(cores.SEAT10_4.depth, mean(cores.SEAT10_4.ages,2) - ...
    std(cores.SEAT10_4.ages,[],2), 'b--')

% Age plots for SEAT2010-5
trace = round(length(radar5.Easting)/2);
age5 = mean(squeeze(radar5.ages(:,trace,:)),2);
std5 = std(squeeze(radar5.ages(:,trace,:)),[],2);
figure
hold on
plot(radar5.depth, age5, 'r', 'LineWidth', 2)
plot(radar5.depth, age5+std5, 'r--')
plot(radar5.depth, age5-std5, 'r--')
plot(cores.SEAT10_5.depth, mean(cores.SEAT10_5.ages,2),...
    'b', 'LineWidth', 2)
plot(cores.SEAT10_5.depth, mean(cores.SEAT10_5.ages,2) + ...
    std(cores.SEAT10_5.ages,[],2), 'b--')
plot(cores.SEAT10_5.depth, mean(cores.SEAT10_5.ages,2) - ...
    std(cores.SEAT10_5.ages,[],2), 'b--')

% Age plots for SEAT2010-6
trace = round(length(radar6.Easting)/2);
age6 = mean(squeeze(radar6.ages(:,trace,:)),2);
std6 = std(squeeze(radar6.ages(:,trace,:)),[],2);
figure
hold on
plot(radar6.depth, age6, 'r', 'LineWidth', 2)
plot(radar6.depth, age6+std6, 'r--')
plot(radar6.depth, age6-std6, 'r--')
plot(cores.SEAT10_6.depth, mean(cores.SEAT10_6.ages,2),...
    'b', 'LineWidth', 2)
plot(cores.SEAT10_6.depth, mean(cores.SEAT10_6.ages,2) + ...
    std(cores.SEAT10_6.ages,[],2), 'b--')
plot(cores.SEAT10_6.depth, mean(cores.SEAT10_6.ages,2) - ...
    std(cores.SEAT10_6.ages,[],2), 'b--')