% Script to help test the sensitivity of of annual SMB estimates to logistic
% regression parameters for all three core sites

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
                %                 data_path = 'C:/Users/durba/Documents/Research/Antarctica/Data/';
                data_path = 'F:/Research/Antarctica/Data/';
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
% Add CReSIS OIB MATLAB reader functions to path
addon_folder = fullfile(addon_path, 'cresis-L1B-matlab-readers/');
addpath(genpath(addon_folder))

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Define number of Monte Carlo simulations to perform
Ndraw = 100;

horz_res = 25;

%%

sites = {'SEAT10_4', 'SEAT10_5', 'SEAT10_6'};

err_up = zeros(1, length(sites));
err_down = zeros(1,length(sites));

for i = 1:length(sites)
    
    % Name of current core site
    site_nm = sites{i};
    
    % Load relevant radar data (previously generated using the above section)
    radar = load(fullfile(data_path, 'IceBridge/manual_layers', site_nm, ...
        strcat('layers_', site_nm, '_opt.mat')));
    
    %Preallocate cell array for auto layer position subscripts
    layers = cell(1, max(radar.groups(:)));
    for j = 1:length(layers)
        
        layers{j} = find(radar.groups==j);
    end
    
    % Remove empty cells from layer array
    layers = layers(~cellfun(@isempty,layers));
    
    % Calculate continuous layer distances for each layer (accounting for
    % lateral size of stacked radar trace bins)
    layers_dist = cellfun(@(x) numel(x)*horz_res, layers);
    
    % Map layer prominence-distance values to the location within the radar
    % matrix of the ith layer
    layer_peaks = zeros(size(radar.peaks));
    for j = 1:length(layers)
        layer_peaks(layers{j}) = radar.peaks(layers{j}).*layers_dist(j);
    end
    
    % Define surface age and the year associated with the first pick of the
    % algorithm
    age_top = radar.collect_date;
    yr_pick1 = ceil(radar.collect_date - 1);
    
    
    %%
    
    
    % Preallocate arrays for layer likelihoods and anges
    ages = zeros([size(radar.data_smooth) Ndraw]);
    likelihoods = zeros(size(radar.data_smooth));
    err_out = [];
    for j = 1:size(layer_peaks, 2)
        
        % Get layer prom-distance values and depths for layers in ith trace
        peaks_i = layer_peaks(:,j);
        peaks_idx = peaks_i>0;
        peaks_i = peaks_i(peaks_idx);
        depths_i = radar.depth(peaks_idx);
        
        % Likelihood of layer representing a year based on a logistic function
        r = -2.4333e-4; % [-3.18e-4 -1.55e-4]
        k = 4.4323;     % [3.25 4.8]
        
        likelihood = 1./(1+exp(r*peaks_i + k));
        likelihoods(peaks_idx,j) = likelihood;
        
        % Assign MC simulation annual layer presence based on layer likelihood
        % values
        yr_idx = zeros(length(depths_i), Ndraw);
        for k = 1:length(depths_i)
            R = rand(Ndraw, 1) <= likelihood(k);
            yr_idx(k,:) = R;
        end
        
        for k = 1:Ndraw
            depths_k = [0; depths_i(logical(yr_idx(:,k)))];
            yrs_k = ([age_top yr_pick1:-1:yr_pick1-length(depths_k)+2])';
            try
                ages(:,j,k) = interp1(depths_k, yrs_k, radar.depth, 'linear', 'extrap');
            catch
                sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
                err_out = [err_out k];
            end
        end
        if ~isempty(err_out)
            ages(:,j,err_out) = repmat(sum(squeeze(ages(:,j,:)), 2)./...
                sum(squeeze(ages(:,j,:))~=0, 2), 1, length(err_out));
        end
        err_out = [];
    end
    
    
    
    % Preallocate arrays for layer likelihoods and anges
    ages_min = zeros([size(radar.data_smooth) Ndraw]);
    likelihoods_min = zeros(size(radar.data_smooth));
    err_out = [];
    for j = 1:size(layer_peaks, 2)
        
        % Get layer prom-distance values and depths for layers in ith trace
        peaks_i = layer_peaks(:,j);
        peaks_idx = peaks_i>0;
        peaks_i = peaks_i(peaks_idx);
        depths_i = radar.depth(peaks_idx);
        
        % Likelihood of layer representing a year based on a logistic function
        r = -3.18e-4; % [-3.18e-4 -1.55e-4]
        k = 3.25;     % [3.25 4.8]
        
        likelihood = 1./(1+exp(r*peaks_i + k));
        likelihoods_min(peaks_idx,j) = likelihood;
        
        % Assign MC simulation annual layer presence based on layer likelihood
        % values
        yr_idx = zeros(length(depths_i), Ndraw);
        for k = 1:length(depths_i)
            R = rand(Ndraw, 1) <= likelihood(k);
            yr_idx(k,:) = R;
        end
        
        for k = 1:Ndraw
            depths_k = [0; depths_i(logical(yr_idx(:,k)))];
            yrs_k = ([age_top yr_pick1:-1:yr_pick1-length(depths_k)+2])';
            try
                ages_min(:,j,k) = interp1(depths_k, yrs_k, radar.depth, 'linear', 'extrap');
            catch
                sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
                err_out = [err_out k];
            end
        end
        if ~isempty(err_out)
            ages_min(:,j,err_out) = repmat(sum(squeeze(ages_min(:,j,:)), 2)./...
                sum(squeeze(ages_min(:,j,:))~=0, 2), 1, length(err_out));
        end
        err_out = [];
    end
    
    
    % Preallocate arrays for layer likelihoods and anges
    ages_max = zeros([size(radar.data_smooth) Ndraw]);
    likelihoods_max = zeros(size(radar.data_smooth));
    err_out = [];
    for j = 1:size(layer_peaks, 2)
        
        % Get layer prom-distance values and depths for layers in ith trace
        peaks_i = layer_peaks(:,j);
        peaks_idx = peaks_i>0;
        peaks_i = peaks_i(peaks_idx);
        depths_i = radar.depth(peaks_idx);
        
        % Likelihood of layer representing a year based on a logistic function
        r = -1.55e-4; % [-3.18e-4 -1.55e-4]
        k = 4.8;     % [3.25 4.8]
        
        likelihood = 1./(1+exp(r*peaks_i + k));
        likelihoods_max(peaks_idx,j) = likelihood;
        
        % Assign MC simulation annual layer presence based on layer likelihood
        % values
        yr_idx = zeros(length(depths_i), Ndraw);
        for k = 1:length(depths_i)
            R = rand(Ndraw, 1) <= likelihood(k);
            yr_idx(k,:) = R;
        end
        
        for k = 1:Ndraw
            depths_k = [0; depths_i(logical(yr_idx(:,k)))];
            yrs_k = ([age_top yr_pick1:-1:yr_pick1-length(depths_k)+2])';
            try
                ages_max(:,j,k) = interp1(depths_k, yrs_k, radar.depth, 'linear', 'extrap');
            catch
                sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
                err_out = [err_out k];
            end
        end
        if ~isempty(err_out)
            ages_max(:,j,err_out) = repmat(sum(squeeze(ages_max(:,j,:)), 2)./...
                sum(squeeze(ages_max(:,j,:))~=0, 2), 1, length(err_out));
        end
        err_out = [];
    end
    
    
    
    
  %%
  
  bias_max = mean(ages_max,3) - mean(ages,3);
  bias_min = mean(ages_min,3) - mean(ages,3);
  
%   figure
%   hold on
%   plot(radar.depth, mean(bias_max,2), 'r')
%   plot(radar.depth, mean(bias_max,2) + std(bias_max,[],2), 'r--')
%   plot(radar.depth, mean(bias_max,2) - std(bias_max,[],2), 'r--')
%   plot(radar.depth, mean(bias_min,2), 'b')
%   plot(radar.depth, mean(bias_min,2) + std(bias_min,[],2), 'b--')
%   plot(radar.depth, mean(bias_min,2) - std(bias_min,[],2), 'b--')
%   hold off
  
  p_max = polyfit(radar.depth, mean(bias_max,2), 1);
  err_up(i) = p_max(1);
  p_min = polyfit(radar.depth, mean(bias_min,2), 1);
  err_down(i) = p_min(1);
  
    
end

err_tot = mean([abs(err_up) abs(err_down)]);
