% Script to optimize logistic regression coefficients for layer-likelihood
% calculations within the PAIPR algorithms


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

% Add PAIPR optimization functions to path
opt_path = fullfile(cd, '..', 'optimize');
addpath(genpath(opt_path))

% Select path to directoryin containing training datasets
path_new = fullfile(data_path, "PAIPR-results/v0.3.0/optim");
% path_new = uigetdir(data_path, ...
%     'Select optimization directory to use');

%% Process raw OIB echograms with PAIPR for use in log regression optimization

% Define number of Monte Carlo simulations to perform
Ndraw = 100;


% Select modeled density subset .mat file to use, and load to workspace
rho_fn = "rho_20111109subset.mat";
rho_path = "/media/durbank/WARP/Research/Antarctica/Data/PAIPR-results/v0.3.0/";
% [rho_fn, rho_path] = uigetfile([data_path '*.mat'], ...
%     "Select modeled density .mat subset");
rho_subset = load(fullfile(rho_path, rho_fn));
rho_subset = rho_subset.rho_subset;


% Define directory containing raw echogram files
core_site = 'SEAT2010_4';
echo_dir = fullfile(path_new, 'raw_data', core_site);


% Get grouping indices of echograms
[files, start_idx, end_idx] = OIB_chunk(echo_dir);






horz_res = 25;


echo_0 = orderfields(import_radar(...
    OIB_convert(fullfile(files(1).folder, files(1).name))));
fields = fieldnames(echo_0);
for i = 1:length(fields)
    echo_0.(fields{i}) = [];
end



for i=1:length(end_idx)
    
    
    echo_i = echo_0;
    for j = start_idx(i):end_idx(i)
        
        OIB_j = OIB_convert(fullfile(files(j).folder, files(j).name));
        struct_j = orderfields(import_radar(OIB_j));
        
        echo_i = cell2struct(cellfun(@horzcat, struct2cell(echo_i), ...
            struct2cell(struct_j), 'uni', 0), fieldnames(echo_i), 1);
    
        
    end
    echo_i.dist = pathdistps(echo_i.lat, echo_i.lon);
    
    
    
    [radar] = radar_stack(echo_i, horz_res);
    
    
    % Load modeled depth-density data from stats model output at specified
    % Easting/Northing coordinates
    [rho_data] = load_rho(rho_subset, radar.Easting, radar.Northing);
    
    % Convert to depth
    [radar] = radar_depth(radar, rho_data);
    
    % Calculate radar age-depth profile distributions (includes processing
    %signal-noise, radon transforms, layer tracing, likelihood assignments,
    % and age calculations)
    [radar] = calc_layers(radar, 'stream');
    
    
    
    
    
    
    
    
    man_layers = load(fullfile(path_new, 'man_layers', core_site));
    man_layers = man_layers.man_layers;
    depth_max = ceil(max(cellfun(@(x) max(x(:,2)), man_layers)));
    dist_max = ceil(max(cellfun(@(x) max(x(:,1)), man_layers)));
    
    man_idx = cellfun(@(x) sub2ind([depth_max, dist_max], ...
        round(x(:,2)),x(:,1)), man_layers, 'UniformOutput', false);
    man_grid = zeros([depth_max dist_max]);
    for k = 1:length(man_idx)
        man_grid(man_idx{k}) = 1;
    end
    
    
    max_idx = min([depth_max length(radar.depth)]);
    
    radar.depth = radar.depth(1:max_idx);
    radar.data_smooth = radar.data_smooth(1:max_idx,:);
    radar.IM_grad = radar.IM_grad(1:max_idx,:);
    radar.IM_QC = radar.IM_QC(1:max_idx,:);
    radar.DB = radar.DB(1:max_idx,:);
    radar.man_grid = man_grid(1:max_idx,1:length(radar.dist));
    
    
    save(fullfile(path_new, "interim_data", strcat(core_site, ".mat")), ...
        '-struct', 'radar')
    
end




%%
clearvars -except data_path path_new 


core_site = "SEAT2010_4";

radar = load(fullfile(path_new, "interim_data", core_site));

DB = radar.DB;
depth = radar.depth;

r_params = zeros(1, size(DB,2));
k_params = zeros(1, size(DB,2));
SSE = zeros(1, size(DB,2));

parfor i = 1:length(r_params)
    
    [r_params(i), k_params(i), SSE(i)] = opt_param(radar.man_grid(:,i), ...
        DB(:,i), depth, 3, -8);
    
end


%%% Save the parameter distributions

output_path = fullfile(path_new, strcat("params_", core_site, ".mat"));
save(output_path, 'r_params', 'k_params', 'SSE')