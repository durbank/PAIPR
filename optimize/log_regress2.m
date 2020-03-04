% Script to estimate the logistic regression coefficients for calculating
% layer likelihood scores

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'G:/Research/Antarctica/Data/';
        addon_path = 'N:/MATLAB/Add-ons/';
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile(addon_path, 'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add export_fig to path
addon_struct = dir(fullfile(addon_path, 'altmany-export_fig*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add CReSIS OIB MATLAB reader functions to path
addon_struct = dir(fullfile(addon_path, 'cresis-L1B-matlab-readers*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add PAIPR-core functions to path
PAIPR_path = fullfile(cd,'..', 'PAIPR-core');
addpath(genpath(PAIPR_path))


% path_new = uigetdir(data_path, ...
%     'Select optimization directory to use');
path_new = fullfile(data_path, "PAIPR-results/tmp/optim");

%% Process raw OIB echograms with PAIPR for use in log regression optimization

% % Define number of Monte Carlo simulations to perform
% Ndraw = 100;
% 
% 
% radar_old = load(fullfile(path_new, "interim_data/interim_PIG.mat"));
% radar = radar_old;
% vert_res = round(mean(diff(radar.depth)),2);
% horz_res = round(mean(diff(radar.dist)));
% 
% 
% % Iterative radon transforms
% [IM_gradients] = radar_gradient(radar, vert_res, horz_res);
% 
% % Find radar peaks in echogram
% [peaks_raw, peak_width] = radar_peaks(radar, vert_res);
% 
% 
% 
% [peaks, group_num, layers] = radar_trace(peaks_raw, peak_width, ...
%     IM_gradients, vert_res, horz_res);
% 
% % Find max depth with full layer coverage
% max_depth = zeros(1,size(group_num,2));
% for i=1:length(max_depth)
%     max_depth(i) = find(group_num(:,i), 1, 'last');
% end
% cut_idx = min(max_depth);
% 
% [rows,cols] = cellfun(@(x) ind2sub(size(peaks),x), layers, ...
%     'UniformOutput', false);
% [layers_new] = cellfun(@(x,y) sub2ind([cut_idx length(radar.dist)], ...
%     x(x<=cut_idx), y(x<=cut_idx)), rows, cols, 'UniformOutput', false);
% layers_new = layers_new(~cellfun(@isempty, layers_new));
% 
% 
% radar.depth = radar.depth(1:cut_idx);
% radar.data_smooth = radar.data_smooth(1:cut_idx,:);
% radar.IM_grad = IM_gradients(1:cut_idx,:);
% radar.peaks = peaks(1:cut_idx,:);
% radar.layers = layers_new;
% radar.groups = group_num(1:cut_idx,:);
% 
% 
% [man_new] = cellfun(@(x) x(x(:,2)<=cut_idx,:), radar.man_layers, ...
%     'UniformOutput', false);
% man_new = man_new(~cellfun(@isempty, man_new));
% 
% radar.man_layers = man_new;
% 
% %%
% % Calculate continuous layer distances for each layer (accounting for 
% % lateral size of stacked radar trace bins)
% layers_dist = cellfun(@(x) numel(x)*horz_res, radar.layers);
% 
% % Find the mean peak prominence (used to scale the prominence-distance
% % results)
% peak_w = 1/std(radar.data_smooth(:));
% dist_w = 1/(size(radar.data_smooth,2)*horz_res);
% 
% % Map layer prominence-distance values to the location within the radar
% % matrix of the ith layer
% layer_DB = zeros(size(radar.peaks));
% for i = 1:length(radar.layers)
%     layer_DB(radar.layers{i}) = peak_w*dist_w*...
%         radar.peaks(radar.layers{i}).*layers_dist(i);
% end
% radar.DB = layer_DB;
% 
% %% Diagnostic plots
% 
% % figure
% % imagesc(radar.data_smooth)
% % hold on
% % for i=1:length(radar.layers)
% %     [r,c] = cellfun(@(x) ind2sub(size(radar.groups), x), ...
% %         radar.layers, 'UniformOutput', false);
% %     cellfun(@(x,y) plot(x, y, 'LineWidth', 2), c, r);
% % end
% % cellfun(@(x) plot(x(:,1), x(:,2), 'k--'), radar.man_layers)
% 
% 
% 
% % Save processed radar structure for future use
% save(fullfile(path_new, "PAIPR_PIG.mat"), '-struct', 'radar')


%% Run `draw_manual.m` to manually trace layers and add to radar structure


%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Load PAIPR radar data
radar_file = fullfile(path_new, "PAIPR_PIG.mat");
radar = load(radar_file);

%%

man_idx = cellfun(@(x) sub2ind(size(radar.data_smooth), ...
    round(x(:,2)),x(:,1)), radar.man_layers, 'UniformOutput', false);
man_grid = zeros(size(radar.data_smooth));
for i = 1:length(man_idx)
    man_grid(man_idx{i}) = 1;
end


DB = radar.DB;
depth = radar.depth;

r_params = zeros(1, size(DB,2));
k_params = zeros(1, size(DB,2));
SSE = zeros(1, size(DB,2));

parfor i = 1:length(r_params)
    
    [r_params(i), k_params(i), SSE(i)] = opt_param(man_grid(:,i), ...
        DB(:,i), depth, 3, -8);
    
end
%% Save the parameter distributions

% output_path = fullfile(path_new, "params_PIG.mat");
% save(output_path, 'r_params', 'k_params', 'SSE')
