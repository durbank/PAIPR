% Script to compare the effect of differences in the assigned value of P_50

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'work';
        %         computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
                
            case 'laptop'
                data_path = 'C:/Users/durba/Documents/Research/Antarctica/Data/';
                addon_path = 'C:/Users/durba/Documents/MATLAB/Addons/';
        end
        
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = fullfile(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))
% Add export_fig to path
addon_folder = fullfile(addon_path, 'altmany-export_fig-cafc7c5/');
addpath(genpath(addon_folder))

%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

Ndraw = 100;

% Load previously processed full OIB radar transect
radar_full = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/OIB_SEAT10_4to10_6_full.mat'));

%% Calculate ages based on different P_50 values (modified from radar_age.m)
tic
% Define horizontal resolution (25 m)
horz_res = 25;

% Calculate continuous layer distances for each layer (accounting for 
% lateral size of stacked radar trace bins)
layers_dist = cellfun(@(x) numel(x)*horz_res, radar_full.layers);

% Map layer prominence-distance values to the location within the radar
% matrix of the ith layer
layer_peaks = zeros(size(radar_full.peaks));
for i = 1:length(radar_full.layers)
    layer_peaks(radar_full.layers{i}) = radar_full.peaks(radar_full.layers{i}).*layers_dist(i);
end

% Define surface age and the year associated with the first pick of the 
% algorithm
age_top = radar_full.collect_date;
yr_pick1 = ceil(radar_full.collect_date - 1);

ages_500 = zeros([size(radar_full.data_smooth) Ndraw]);
err_out = [];
P_50_dist = 500;
for i = 1:size(layer_peaks, 2)
    
    % Assign the 50% likelihood point based on median trace prominence and
    % layer length
    P_50 = median(radar_full.peaks(:,i))*min([P_50_dist 0.5*radar_full.dist(end)]);
%     P_50 = median(Proms{i})*mean(cellfun(@length, layers_idx));
    
    % Assign min/max layer likelihoods, and calculate the logistic rate
    % coefficient
    Po = 0.05;
    K = 1;
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = layer_peaks(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar_full.depth(peaks_idx);
    
    % Likelihood of layer representing a year based on a logistic function
    % with rate (r) calculated above
    likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i));
    
    % Assign MC simulation annual layer presence based on layer likelihood
    % values
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end
    
    for j = 1:Ndraw
        depths_j = [0; depths_i(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        try
            ages_500(:,i,j) = interp1(depths_j, yrs_j, radar_full.depth, 'linear', 'extrap');
        catch
            sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
            err_out = [err_out j];
        end
    end
    if ~isempty(err_out)
        ages_500(:,i,err_out) = repmat(sum(squeeze(ages_500(:,i,:)), 2)./...
            sum(squeeze(ages_500(:,i,:))~=0, 2), 1, length(err_out));
    end
    err_out = [];
end
ages_500 = median(ages_500, 3);

%%% Distance for P_50 is 1 km

ages_1km = zeros([size(radar_full.data_smooth) Ndraw]);
err_out = [];
P_50_dist = 1000;
for i = 1:size(layer_peaks, 2)
    
    % Assign the 50% likelihood point based on median trace prominence and
    % layer length
    P_50 = median(radar_full.peaks(:,i))*min([P_50_dist 0.5*radar_full.dist(end)]);
%     P_50 = median(Proms{i})*mean(cellfun(@length, layers_idx));
    
    % Assign min/max layer likelihoods, and calculate the logistic rate
    % coefficient
    Po = 0.05;
    K = 1;
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = layer_peaks(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar_full.depth(peaks_idx);
    
    % Likelihood of layer representing a year based on a logistic function
    % with rate (r) calculated above
    likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i));
    
    % Assign MC simulation annual layer presence based on layer likelihood
    % values
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end
    
    for j = 1:Ndraw
        depths_j = [0; depths_i(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        try
            ages_1km(:,i,j) = interp1(depths_j, yrs_j, radar_full.depth, 'linear', 'extrap');
        catch
            sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
            err_out = [err_out j];
        end
    end
    if ~isempty(err_out)
        ages_1km(:,i,err_out) = repmat(sum(squeeze(ages_1km(:,i,:)), 2)./...
            sum(squeeze(ages_1km(:,i,:))~=0, 2), 1, length(err_out));
    end
    err_out = [];
end
ages_1km = median(ages_1km, 3);

%%% Distance for P_50 is 10 km

ages_10km = zeros([size(radar_full.data_smooth) Ndraw]);
err_out = [];
P_50_dist = 10000;
for i = 1:size(layer_peaks, 2)
    
    % Assign the 50% likelihood point based on median trace prominence and
    % layer length
    P_50 = median(radar_full.peaks(:,i))*min([P_50_dist 0.5*radar_full.dist(end)]);
%     P_50 = median(Proms{i})*mean(cellfun(@length, layers_idx));
    
    % Assign min/max layer likelihoods, and calculate the logistic rate
    % coefficient
    Po = 0.05;
    K = 1;
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = layer_peaks(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar_full.depth(peaks_idx);
    
    % Likelihood of layer representing a year based on a logistic function
    % with rate (r) calculated above
    likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i));
    
    % Assign MC simulation annual layer presence based on layer likelihood
    % values
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end
    
    for j = 1:Ndraw
        depths_j = [0; depths_i(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        try
            ages_10km(:,i,j) = interp1(depths_j, yrs_j, radar_full.depth, 'linear', 'extrap');
        catch
            sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
            err_out = [err_out j];
        end
    end
    if ~isempty(err_out)
        ages_10km(:,i,err_out) = repmat(sum(squeeze(ages_10km(:,i,:)), 2)./...
            sum(squeeze(ages_10km(:,i,:))~=0, 2), 1, length(err_out));
    end
    err_out = [];
end
ages_10km = median(ages_10km, 3);

%%% Distance for P_50 is 10 km

ages_25km = zeros([size(radar_full.data_smooth) Ndraw]);
err_out = [];
P_50_dist = 25000;
for i = 1:size(layer_peaks, 2)
    
    % Assign the 50% likelihood point based on median trace prominence and
    % layer length
    P_50 = median(radar_full.peaks(:,i))*min([P_50_dist 0.5*radar_full.dist(end)]);
%     P_50 = median(Proms{i})*mean(cellfun(@length, layers_idx));
    
    % Assign min/max layer likelihoods, and calculate the logistic rate
    % coefficient
    Po = 0.05;
    K = 1;
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = layer_peaks(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar_full.depth(peaks_idx);
    
    % Likelihood of layer representing a year based on a logistic function
    % with rate (r) calculated above
    likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i));
    
    % Assign MC simulation annual layer presence based on layer likelihood
    % values
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end
    
    for j = 1:Ndraw
        depths_j = [0; depths_i(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        try
            ages_25km(:,i,j) = interp1(depths_j, yrs_j, radar_full.depth, 'linear', 'extrap');
        catch
            sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
            err_out = [err_out j];
        end
    end
    if ~isempty(err_out)
        ages_25km(:,i,err_out) = repmat(sum(squeeze(ages_25km(:,i,:)), 2)./...
            sum(squeeze(ages_25km(:,i,:))~=0, 2), 1, length(err_out));
    end
    err_out = [];
end
ages_25km = median(ages_25km, 3);
toc

%% 
clearvars -except radar_full ages* cores
ages_5km = median(radar_full.ages, 3);

figure
hold on
h1 = plot(ages_5km(end,:), 'b', 'LineWidth', 2);
% plot(median(squeeze(radar_full.ages(end,:,:)), 2) + ...
%     2*std(squeeze(radar_full.ages(end,:,:)), [], 2), 'b--', 'LineWidth', 0.5)
% plot(median(squeeze(radar_full.ages(end,:,:)), 2) - ...
%     2*std(squeeze(radar_full.ages(end,:,:)), [], 2), 'b--', 'LineWidth', 0.5)
h2 = plot(ages_500(end,:), 'r', 'LineWidth', 2);
% plot(median(squeeze(ages_500(end,:,:)), 2) + ...
%     2*std(squeeze(ages_500(end,:,:)), [], 2), 'r--', 'LineWidth', 0.5)
% plot(median(squeeze(ages_500(end,:,:)), 2) - ...
%     2*std(squeeze(ages_500(end,:,:)), [], 2), 'r--', 'LineWidth', 0.5)
h3 = plot(ages_1km(end,:), 'm', 'LineWidth', 2);
% plot(median(squeeze(ages_1km(end,:,:)), 2) + ...
%     2*std(squeeze(ages_1km(end,:,:)), [], 2), 'm--', 'LineWidth', 0.5)
% plot(median(squeeze(ages_1km(end,:,:)), 2) - ...
%     2*std(squeeze(ages_1km(end,:,:)), [], 2), 'm--', 'LineWidth', 0.5)
h4 = plot(ages_10km(end,:), 'c', 'LineWidth', 2);
% plot(median(squeeze(ages_10km(end,:,:)), 2) + ...
%     2*std(squeeze(ages_10km(end,:,:)), [], 2), 'c--', 'LineWidth', 0.5)
% plot(median(squeeze(ages_10km(end,:,:)), 2) - ...
%     2*std(squeeze(ages_10km(end,:,:)), [], 2), 'c--', 'LineWidth', 0.5)
h5 = plot(ages_25km(end,:), 'g', 'LineWidth', 2);