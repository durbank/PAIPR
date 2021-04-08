% Temporary file for testing new features for process_PAIPR.m, prior to
% migrating those features to the main file

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/Data/';
        addon_path = 'N:/MATLAB/Add-Ons';
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
addon_struct = dir(fullfile(addon_path, 'cresis-L1B-matlab-reader*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Import PAIPR functions
addpath(genpath("../PAIPR-core"))

%%

% Define number of Monte Carlo simulations to perform
Ndraw = 100;

% Define path to the directory containing radar data (relative to the
% 'data' directory path)
radar_dir = uigetdir(data_path, ...
    'Select folder containing .mat radar data files');

% Define path to output directory of where to store results
output_dir = uigetdir(data_path, ...
    "Select folder to store output .mat results files");
gamma_dir = uigetdir(data_path, ...
    "Select folder to store output gamma-fitted .csv results");

%%

% Select modeled density subset .mat file to use, and load to workspace
[rho_fn, rho_path] = uigetfile([data_path '*.mat'], ...
    "Select modeled density .mat subset");
rho_subset = load(fullfile(rho_path, rho_fn));
rho_subset = rho_subset.rho_subset;

%%

% Get grouping indices of echograms
[files, start_idx, end_idx] = OIB_chunk(radar_dir);


% Define overlap distance and the final horizontal resolution of the output
% data
overlap = 10000;
horz_res = 25;


echo_0 = orderfields(import_radar(...
    OIB_convert(fullfile(files(1).folder, files(1).name))));
fields = fieldnames(echo_0);
for i = 1:length(fields)
    echo_0.(fields{i}) = [];
end



parfor i=1:length(end_idx)
    
    
    echo_i = echo_0;
    for j = start_idx(i):end_idx(i)
        
        OIB_j = OIB_convert(fullfile(files(j).folder, files(j).name));
        struct_j = orderfields(import_radar(OIB_j));
        
        echo_i = cell2struct(cellfun(@horzcat, struct2cell(echo_i), ...
            struct2cell(struct_j), 'uni', 0), fieldnames(echo_i), 1);
        
        
    end
    echo_i.dist = pathdistps(echo_i.lat, echo_i.lon);
    
    
    
    [radar_tmp] = radar_stack(echo_i, horz_res);
    
    
    % Load modeled depth-density data from stats model output at specified
    % Easting/Northing coordinates
    [rho_data] = load_rho(rho_subset, radar_tmp.Easting, radar_tmp.Northing);
    
    % Convert to depth
    [radar_tmp] = radar_depth(radar_tmp, rho_data);
    
    % Calculate radar age-depth profile distributions (includes processing
    %signal-noise, radon transforms, layer tracing, likelihood assignments,
    % and age calculations)
    [radar_tmp] = calc_layers(radar_tmp, 'stream');
    
    % Calculate age-depth profiles
    r = -6.723;
    k = 3.0;
    [radar_tmp] = radar_age(radar_tmp, r, k, Ndraw);
    
    % Calculate radar annual SMB
    [radar_tmp] = calc_SWE(radar_tmp, rho_data, Ndraw);
    
    % Perform QC check on echogram image
    [QC_med, QC_val, QC_flag, depth_idx] = QC_check(...
    radar_tmp.IM_QC, 0.50, 0.10);
    
    % Clip radar data structure variables based on the desired radargram
    % overlap, and combine desired clipped variables into new data
    % structure
    clip = round(0.5*overlap/horz_res);
    radar = struct('collect_time',radar_tmp.collect_time(clip:end-clip),...
        'Easting', radar_tmp.Easting(clip:end-clip),...
        'Northing', radar_tmp.Northing(clip:end-clip), ...
        'dist', radar_tmp.dist(clip:end-clip), 'depth', radar_tmp.depth, ...
        'data_smooth', radar_tmp.data_smooth(:,clip:end-clip),...
        'IM_grad', radar_tmp.IM_grad(:,clip:end-clip), ...
        'IM_QC', radar_tmp.IM_QC(:,clip:end-clip), 'QC_med', QC_med, ...
        'QC_val', QC_val, 'QC_flag', QC_flag', 'QC_depth_idx', depth_idx,...
        'DB', radar_tmp.DB(:,clip:end-clip), ...
        'likelihood', radar_tmp.likelihood(:,clip:end-clip), ...
        'ages', radar_tmp.ages(:,clip:end-clip,:));
    
    
    if isfield(radar_tmp, 'groups')
        radar.groups = radar_tmp.groups(:,clip:end-clip);
        
        % Find layer member indices based on new clipped record
        layers = cell(1, max(radar.groups(:)));
        for j = 1:length(layers)
            layers{j} = find(radar.groups == j);
        end
        radar.layers = layers(~cellfun(@isempty, layers));
        
    end
    
    % Check for existing elevation data, and add to output struct if
    % available
    if isfield(radar_tmp, 'elev')
        radar.elev = radar_tmp.elev(clip:end-clip);
    else
        radar.elev = nan(1, length(radar.Easting));
    end
    
    % Redefine radargram distances based on clipped data
    radar.dist = radar.dist - radar.dist(1);
    
    % Clip SMB cell array data and add to output struct
    radar.SMB_yr =  radar_tmp.SMB_yr(clip:end-clip);
    radar.SMB = radar_tmp.SMB(clip:end-clip);
    
    % Generate file names and paths under which to save data
    filename = sprintf('%s%d','radar_out',i);
    mat_output = fullfile(output_dir, strcat(filename, '.mat'));
    csv_output = fullfile(gamma_dir, strcat(filename, '.csv'));
    
    % Save output structures to disk
    [save_success1] = parsave_all(radar, mat_output, csv_output);
    
end