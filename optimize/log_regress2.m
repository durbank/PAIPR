% Script to estimate the logistic regression coefficients for calculating
% layer likelihood scores

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
% Add export_fig to path
addon_struct = dir(fullfile(addon_path, 'altmany-export_fig*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add CReSIS OIB MATLAB reader functions to path
addon_struct = dir(fullfile(addon_path, 'cresis-L1B-matlab-readers*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add PAIPR-core functions to path
PAIPR_path = fullfile(cd,'..', 'PAIPR-core');
addpath(genpath(PAIPR_path))


path_new = uigetdir(data_path, ...
    'Select optimization directory to use');

%% Process raw OIB echograms with PAIPR for use in log regression optimization

% Define number of Monte Carlo simulations to perform
Ndraw = 100;


radar_old = load(fullfile(path_new, "PAIPR_SEAT6.mat"));
radar = radar_old;
vert_res = round(mean(diff(radar.depth)),2);
horz_res = round(mean(diff(radar.dist)));



% Iterative radon transforms
[IM_gradients] = radar_gradient(radar, vert_res, horz_res);

% Find radar peaks in echogram
[peaks_raw, peak_width] = radar_peaks(radar, vert_res);



[peaks, group_num, layers] = radar_trace(peaks_raw, peak_width, ...
    IM_gradients, vert_res, horz_res);

% Find max depth with full layer coverage
max_depth = zeros(1,size(group_num,2));
for i=1:length(max_depth)
    max_depth(i) = find(group_num(:,i), 1, 'last');
end
cut_idx = min(max_depth);

[rows,cols] = cellfun(@(x) ind2sub(size(peaks),x), layers, ...
    'UniformOutput', false);
[layers_new] = cellfun(@(x,y) sub2ind([cut_idx length(radar.dist)], ...
    x(x<=cut_idx), y(x<=cut_idx)), rows, cols, 'UniformOutput', false);
layers_new = layers_new(~cellfun(@isempty, layers_new));


radar.depth = radar.depth(1:cut_idx);
radar.data_smooth = radar.data_smooth(1:cut_idx,:);
radar.IM_grad = IM_gradients(1:cut_idx,:);
radar.peaks = peaks(1:cut_idx,:);
radar.layers = layers_new;
radar.groups = group_num(1:cut_idx,:);




[man_new] = cellfun(@(x) x(x(:,2)<=cut_idx,:), radar.man_layers, ...
    'UniformOutput', false);
man_new = man_new(~cellfun(@isempty, man_new));

radar.man_layers = man_new;

%%
% Save processed radar structure for future use
% save(fullfile(path_new, "PAIPR_PIG.mat"), '-struct', 'radar')


%% Run `draw_manual.m` to manually trace layers and add to radar structure


%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Load PAIPR radar data
radar_file = fullfile(path_new, "PAIPR_SEAT6.mat");
radar = load(radar_file);

%%

man_idx = cellfun(@(x) sub2ind(size(radar.data_smooth), ...
    round(x(:,2)),x(:,1)), radar.man_layers, 'UniformOutput', false);
man_grid = zeros(size(radar.data_smooth));
for i = 1:length(man_idx)
    man_grid(man_idx{i}) = 1;
end

%%

% cut_idx = 1100;
% DB = radar.DB(1:cut_idx,:);
% man_cut = man_grid(1:cut_idx,:);
% depth = radar.depth(1:cut_idx);
DB = radar.DB;
depth = radar.depth;

r_params = zeros(1, size(DB,2));
k_params = zeros(1, size(DB,2));
SSE = zeros(1, size(DB,2));

parfor i = 1:length(r_params)
    
    [r_params(i), k_params(i), SSE(i)] = opt_param(man_grid(:,i), ...
        DB(:,i), depth, 3, -8);
    
end


%%

% output_path = fullfile(path_new, "params_PIG.mat");
% save(output_path, 'r_params', 'k_params', 'SSE')
