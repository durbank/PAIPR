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

% Combine all input echograms and decompose into echogram subdomains of
% desired length
radar_dir = fullfile(path_new, 'raw_input/');
radar_ALL = echo_format(radar_dir);

% Select modeled density subset .mat file to use, and load to workspace
rho_subset = load(fullfile(path_new, 'rho_20111109subset'));
rho_subset = rho_subset.rho_subset;

% Define overlap distance and the final horizontal resolution of the output
% data
overlap = 10000;
horz_res = 25;

radar_tmp = struct();
names = {'SEAT10_6', 'SEAT10_5', 'SEAT10_4'};

for i = 1:length(radar_ALL)
    
    % Find the mean response with depth in the radar data attributes across
    % a given horizontal resolution (in meters)
    [radar] = radar_stack(radar_ALL(i).segment, horz_res);
    
    % Load modeled depth-density data from stats model output at specified
    % Easting/Northing coordinates
    [rho_data] = load_rho(rho_subset, radar.Easting, radar.Northing);
    
    % Convert to depth
    [radar] = radar_depth(radar, rho_data);
    
    %%
    % Calculate radar age-depth profile distributions (includes processing
    %signal-noise, radon transforms, layer tracing, likelihood assignments,
    % and age calculations)
    
    % Stationarize the radar response by differencing traces with a smoothing
    % spline
    s = zeros(size(radar.data_stack));
    for j = 1:size(s, 2)
        s(:,j) = csaps(radar.depth(:,j), radar.data_stack(:,j), ...
            0.95, radar.depth(:,j));
    end
    radar_stat = radar.data_stack - s;
    
    % Define the vertical resolution of the core data and horizontal resolution
    % of the radar data
    vert_res = 0.02;
    horz_res = round(mean(diff(radar.dist)));
    
    % Define the cutoff depth for radar traces and find index of crossover
    % depth
    cutoff = 30;
    depth_bott = floor(min([min(radar.depth(end,:)) cutoff]));
    
    % Trim radar traces to cutoff depth and interpolate data to vertical scale
    % of the firn cores
    radar_interp = zeros(depth_bott/vert_res+1, size(radar.data_stack, 2));
    for j = 1:size(radar.data_stack, 2)
        depth_interp = (0:vert_res:radar.depth(end,j));
        radar_i = interp1(radar.depth(:,j), radar_stat(:,j), ...
            depth_interp, 'pchip');
        radar_interp(:,j) = radar_i(1:size(radar_interp, 1));
    end
    
    % Assign structure output depth to interpolated depths
    radar.depth = (0:vert_res:depth_bott)';
    
    % Smooth the laterally averaged radar traces with depth based on a 3rd
    % order Savitzky-Golay filter with a window of 9 frames (~18 cm)
    radar.data_smooth = sgolayfilt(radar_interp, 3, 9);
    
    
    
    
    % Iterative radon transforms
    [IM_gradients] = radar_gradient(radar, vert_res, horz_res);
    
    %
    [stream_val, XY_stream] = stream_sum(radar, IM_gradients);
    
    
    %
    [~, peak_idx, ~, peak_prom] = findpeaks(stream_val, ...
        'MinPeakProminence', std(stream_val)/10);
    
    
    % Find the mean peak prominence (used to scale the prominence-distance
    % results)
    peak_w = 1/std(radar.data_smooth(:));
    dist_w = 1/(size(radar.data_smooth,2)*horz_res);
    % peak_w = 1;
    % dist_w = 1;
    
    
    layer_DB = zeros(size(radar.data_smooth));
    for j=1:length(peak_idx)
        
        layer_DB(XY_stream{peak_idx(j)}) = peak_w*dist_w*peak_prom(j);
    end
    
    %%
    % Clip radar data structure variables based on the desired radargram
    % overlap, and combine desired clipped variables into new data
    % structure
    radar_tmp.(names{i}) = struct(...
        'collect_time',radar.collect_time, ...
        'Easting', radar.Easting,...
        'Northing', radar.Northing, 'dist', radar.dist, ...
        'depth', radar.depth, 'data_smooth', radar.data_smooth,...
        'DB', layer_DB);
end
radar = radar_tmp;
SEAT4 = radar.SEAT10_4;
SEAT5 = radar.SEAT10_5;
SEAT6 = radar.SEAT10_6;


%
% % Save processed radar structure for future use
% output_path = ;
% output_name = ;
% save(fullfile(output_path, output_name), '-struct', 'radar')


%% Run `draw_manual.m` to manually trace layers and add to radar structure


%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Load PAIPR radar data
radar_file = fullfile(path_new, "PAIPR_SEAT4.mat");
radar = load(radar_file);

%% Calculate distance brightnesses

% Weighting coefficients for mean brightness and echogram distance
peak_w = 1/std(radar.data_smooth(:));
dist_w = 1/radar.dist(end);
% peak_w = 1;
% dist_w = 1;

%%

man_idx = cellfun(@(x) sub2ind(size(radar.data_smooth), ...
    round(x(:,2)),x(:,1)), radar.man_layers, 'UniformOutput', false);
man_grid = zeros(size(radar.data_smooth));
for i = 1:length(man_idx)
    man_grid(man_idx{i}) = 1;
end

%%

cut_idx = 1100;
DB = radar.DB(1:cut_idx,:);
man_cut = man_grid(1:cut_idx,:);
depth = radar.depth(1:cut_idx);

r_params = zeros(1, size(DB,2));
k_params = zeros(1, size(DB,2));
SSE = zeros(1, size(DB,2));

parfor i = 1:length(r_params)
    
    [r_params(i), k_params(i), SSE(i)] = opt_param(man_cut(:,i), ...
        DB(:,i), depth, 3, -8);
    
end


%%

output_path = fullfile(input_dir, "params_output.mat");
save(output_path, 'r_params', 'k_params', 'SSE')
