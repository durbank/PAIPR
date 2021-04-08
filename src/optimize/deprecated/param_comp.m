% Script to compare the results of logistic regression parameter estimation
% between the different validation sites

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/Data/';
        addon_path = 'N:/MATLAB/Add-ons/';
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

% input_dir = uigetdir(data_path, "Select optimization directory");
input_dir = fullfile(data_path, "PAIPR-results/v0.3.0/optim/");

SEAT_4 = load(fullfile(input_dir, "params_SEAT2010_4.mat"));
SEAT_5 = load(fullfile(input_dir, "params_SEAT2010_5.mat"));
SEAT_6 = load(fullfile(input_dir, "params_SEAT2010_6.mat"));
PIG = load(fullfile(input_dir, "params_PIG2010.mat"));

%%
figure
title("Rate parameter (r)")
hold on
ksdensity(SEAT_4.r_params(SEAT_4.SSE<25))
ksdensity(SEAT_5.r_params(SEAT_5.SSE<25))
ksdensity(SEAT_6.r_params(SEAT_6.SSE<25))
ksdensity(PIG.r_params(PIG.SSE<25))
legend("2010-4", "2010-5", "2010-6", "PIG")
hold off

figure
title("Rate parameter (r) SSE")
hold on
ksdensity(SEAT_4.SSE(SEAT_4.SSE<25))
ksdensity(SEAT_5.SSE(SEAT_5.SSE<25))
ksdensity(SEAT_6.SSE(SEAT_6.SSE<25))
ksdensity(PIG.SSE(PIG.SSE<25))
legend("2010-4", "2010-5", "2010-6", "PIG")
hold off

w4 = (1./SEAT_4.SSE)./sum(1./SEAT_4.SSE);
w5 = (1./SEAT_5.SSE)./sum(1./SEAT_5.SSE);
w6 = (1./SEAT_6.SSE)./sum(1./SEAT_6.SSE);
wPIG = (1./PIG.SSE)./sum(1./PIG.SSE);

sprintf("Mean error-weighted r parameters: %0.2f | %0.2f | %0.2f | %0.2f", ...
    sum(w4.*SEAT_4.r_params), sum(w5.*SEAT_5.r_params), ...
    sum(w6.*SEAT_6.r_params), sum(wPIG.*PIG.r_params))

r_all = [SEAT_4.r_params SEAT_5.r_params SEAT_6.r_params PIG.r_params];
SSE_all = [SEAT_4.SSE SEAT_5.SSE SEAT_6.SSE PIG.SSE];
w_all = (1./SSE_all)./sum(1./SSE_all);

r_hat = sum(((1./SSE_all)./sum(1./SSE_all)).*r_all);
r_std = std(r_all, w_all);
sprintf("Final error-weighted estimate of r: %0.3f", r_hat)
sprintf("Finall error-weighted standard deviation of r: %0.3f", r_std)





%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);


[f4, xi4] = ksdensity(SEAT_4.r_params(SEAT_4.r_params < ...
    quantile(SEAT_4.r_params, 0.99) & SEAT_4.r_params > ...
    quantile(SEAT_4.r_params, 0.01)));
[~, f_idx] = max(f4);
r4 = xi4(f_idx);

[f5, xi5] = ksdensity(SEAT_5.r_params(SEAT_5.r_params < ...
    quantile(SEAT_5.r_params, 0.99) & SEAT_5.r_params > ...
    quantile(SEAT_5.r_params, 0.01)));
[~, f_idx] = max(f5);
r5 = xi5(f_idx);

[f6, xi6] = ksdensity(SEAT_6.r_params(SEAT_6.r_params < ...
    quantile(SEAT_6.r_params, 0.99) & SEAT_6.r_params > ...
    quantile(SEAT_6.r_params, 0.01)));
[~, f_idx] = max(f6);
r6 = xi6(f_idx);


figure
hold on
plot(xi4, f4)
plot(xi5, f5)
plot(xi6, f6)


r = mean([r4 r5 r6]);
k = mean([median(SEAT_4.k_params), median(SEAT_5.k_params), ...
    median(SEAT_6.k_params)]);

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


X = 0:max(radar4.DB(:));
% X = 0:0.01:3.5;
Y = 1./(1+exp(r*X + k));
Y_ub = 1./(1+exp((r+2.75e-4)*X + k));
Y_lb = 1./(1+exp((r-2.75e-4)*X + k));

figure
hold on
% plot(radar4.DB(radar4.DB>0), radar4.likelihood(radar4.DB>0), 'o')
% plot(radar5.DB(radar5.DB>0), radar5.likelihood(radar5.DB>0), 'o')
% plot(radar6.DB(radar6.DB>0), radar6.likelihood(radar6.DB>0), 'o')
plot(X,Y, 'b')
plot(X,Y_lb, 'r--')
plot(X,Y_ub, 'r--')


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