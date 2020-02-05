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


path_new = fullfile(data_path, 'PAIPR-results/WAIS_test/optim');

%% Process raw OIB echograms with PAIPR for use in log regression optimization

% Define number of Monte Carlo simulations to perform
Ndraw = 100;

radar_old = load(fullfile(path_new, "Nov25_data/PAIPR_SEAT6.mat"));
vert_res = round(mean(diff(radar_old.depth)),2);
horz_res = round(mean(diff(radar_old.dist)));


% Iterative radon transforms
[IM_gradients] = radar_gradient(radar_old, vert_res, horz_res);

% 
% xstart_idx = 1:300:length(radar.dist);
xstart_idx = round(length(radar_old.dist)/2);
[stream_val, XY_streams] = stream_sum(radar_old, IM_gradients, xstart_idx);

% 
[~, peak_idx, ~, peak_prom] = findpeaks(stream_val, ...
        'MinPeakProminence', std(stream_val)/10);

% Find the mean peak prominence (used to scale the prominence-distance
% results)
peak_w = 1/std(radar_old.data_smooth(:));
dist_w = 1/(size(radar_old.data_smooth,2)*horz_res);


layer_DB = zeros(size(radar_old.data_smooth));
for i=1:length(peak_idx)
    
    layer_DB(XY_streams{peak_idx(i)}) = peak_w*dist_w*peak_prom(i);
end


max_depth = zeros(1,size(layer_DB,2));
for i=1:length(max_depth)
    
    max_depth(i) = find(layer_DB(:,i), 1, 'last');
end

% Clip depth-related variables to final cutoff depth
cut_idx = min(max_depth);



radar = struct('collect_time', radar_old.collect_time, 'Easting', ...
    radar_old.Easting, 'Northing', radar_old.Northing, ...
    'dist', radar_old.dist, 'depth', radar_old.depth(1:cut_idx), ...
    'data_smooth', radar_old.data_smooth(1:cut_idx,:), ...
    'IM_grad', IM_gradients(1:cut_idx,:), 'DB', layer_DB(1:cut_idx,:));

if isfield(radar_old, 'elev')
    radar.elev = radar_old.elev;
end

man_layers = cell(1,length(radar_old.man_layers));
for i=1:length(man_layers)
    
    layer_i = radar_old.man_layers{i};
    keep_idx = layer_i(:,2) <= length(radar.depth);
    man_layers{i} = layer_i(keep_idx,:);
end

radar.man_layers = man_layers(~cellfun(@isempty, man_layers));

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
