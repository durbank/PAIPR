% Wrapper function to import and process an arbitrary number of radargram
% files within a single directory

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'G:/Research/Antarctica/Data/';
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

% % % Load core data from file (data used was previously generated using
% % % import_cores.m)
% % core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
% % cores = load(core_file);
% 
% % Select modeled density file to use, and load in as a table
% [rho_fn, rho_path] = uigetfile([data_path '*.csv'], ...
%     "Select modeled density .csv file");
% rho_file = fullfile(rho_path, rho_fn);
% rho_raw = readtable(rho_file);
% 
% % Convert modelled density lat/lon to Easting/Northing
% [rho_E, rho_N] = ll2ps(rho_raw{:,'Lat'}, rho_raw{:,'Lon'});
% T_locs = table(rho_E, rho_N,'VariableNames', {'Easting', 'Northing'});
% rho_full = [T_locs rho_raw(:,3:end)];
% clear rho_raw rho_E rho_N T_locs

%%

% Combine all input echograms and decompose into echogram subdomains of
% desired length
radar_ALL = echo_format(radar_dir);

% Define overlap distance and the final horizontal resolution of the output
% data
overlap = 10000;
horz_res = 25;

% Determine if all data break segments are sufficiently long to process
fdns = fieldnames(radar_ALL);
for i = 1:length(fdns)
    
    if radar_ALL.(fdns{i}).dist(end) < 1.5*overlap
        radar_ALL = rmfield(radar_ALL, fdns{i});
    end
end
fdns = fieldnames(radar_ALL);
tmp_path = fullfile(data_path, "echo_all_tmp.mat");
save(tmp_path, '-struct', 'radar_ALL')
clear radar_ALL

%%

% depth = 30;         % Depth of modeled density data
% depth_sz = 0.02;    % Vertical resolution of modeled data
% raw_res = 25;       % Horizontal resolution of modeled data
% new_res = 1000;     % Desired horizontal resolution of density data
% depth_num = round(depth/depth_sz);  % Num of observations at each point
% 
% % Indices of starting (surface) locations to select
% select_idx = 1:(new_res/raw_res)*(depth_num+1):size(rho_full,1);
% 
% % Preallocate table for density subset data
% rho_subset = table(0, 0, 'VariableNames', {'Easting', 'Northing'});
% rho_subset = repmat(rho_subset, length(select_idx), 1);
% rho_subset.Data = cell(height(rho_subset),1);
% 
% for i = 1:height(rho_subset)
%     
%     % Add subset data to new table from original table
%     rho_subset.Easting(i) = rho_full.Easting(select_idx(i));
%     rho_subset.Northing(i) = rho_full.Northing(select_idx(i));
%     rho_subset.Data{i} = rho_full(...
%         select_idx(i):(select_idx(i)+depth_num),3:5); 
% end
% 
% clear rho_full

%%

% tmp = load(tmp_path);
% fdns = fieldnames(tmp);
% clear tmp


% Parellel for loop to process all decomposed echograms
parfor i = 1:length(fdns)
    
    radar_tmp = load(tmp_path, fdns{i});
    % Find the mean response with depth in the radar data attributes across
    % a given horizontal resolution (in meters)
    [radar_tmp] = radar_stack(radar_tmp.(fdns{i}), horz_res);
    
    % Load modeled depth-density data from stats model output at specified
    % Easting/Northing coordinates
    [rho_data] = load_rho(rho_subset, radar_tmp.Easting, radar_tmp.Northing);
    
    % Convert to depth
    [radar_tmp] = radar_depth(radar_tmp, rho_data);
    
    % Calculate radar age-depth profile distributions (includes processing
    %signal-noise, radon transforms, layer tracing, likelihood assignments,
    % and age calculations)
    r = -4.3491e-4;
    k = 4.600;
    [radar_tmp] = calc_age(radar_tmp, r, k, Ndraw);
%     [radar_tmp] = calc_age2(radar_tmp, -8.305, 3, Ndraw);
    
    % Calculate radar annual SMB
    [radar_tmp] = calc_SWE(radar_tmp, rho_data, Ndraw);
    
    
    %%
    
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
        'DB', radar_tmp.DB(:,clip:end-clip), ...
        'likelihood', radar_tmp.likelihood(:,clip:end-clip), ...
        'ages', radar_tmp.ages(:,clip:end-clip,:));
    
    
%     % Clip radar data structure variables based on the desired radargram
%     % overlap, and combine desired clipped variables into new data
%     % structure
%     clip = round(0.5*overlap/horz_res);
%     radar = struct('collect_time',radar_tmp.collect_time(clip:end-clip),...
%         'Easting', radar_tmp.Easting(clip:end-clip),...
%         'Northing', radar_tmp.Northing(clip:end-clip), ...
%         'dist', radar_tmp.dist(clip:end-clip), 'depth', radar_tmp.depth, ...
%         'data_smooth', radar_tmp.data_smooth(:,clip:end-clip),...
%         'peaks', radar_tmp.peaks(:,clip:end-clip), ...
%         'groups', radar_tmp.groups(:,clip:end-clip),...
%         'likelihood', radar_tmp.likelihood(:,clip:end-clip), ...
%         'ages', radar_tmp.ages(:,clip:end-clip,:));
%     
%     % Find layer member indices based on new clipped record
%     layers = cell(1, max(radar.groups(:)));
%     for j = 1:length(layers)
%         layers{j} = find(radar.groups == j);
%     end
%     radar.layers = layers(~cellfun(@isempty, layers));
    
    % Redefine radargram distances based on clipped data
    radar.dist = radar.dist - radar.dist(1);
    
    % Check for existing elevation data, and add to output struct if
    % available
    if isfield(radar_tmp, 'elev')
        radar.elev = radar_tmp.elev(clip:end-clip);
    else
        radar.elev = nan(1, length(radar.Easting));
    end
    
    % Clip SMB cell array data and add to output struct
    radar.SMB_yr =  radar_tmp.SMB_yr(clip:end-clip);
    radar.SMB = radar_tmp.SMB(clip:end-clip);
    
    % Generate file names and paths under which to save data
    filename = sprintf('%s%d','radar_out',i);
    mat_output = fullfile(output_dir, strcat(filename, '.mat'));
    csv_output = fullfile(gamma_dir, strcat(filename, '.csv'));
    
    % Save output structures to disk
    [save_success1] = parsave_all(radar, mat_output, csv_output);
    
    %% Output overlapping data for comparisons of neighboring echograms
    
%     % Check for existence of directory for clipped results in output_dir,
%     % and create it if it does not exist
%     clip_dir = "result_clips";
%     if ~exist(fullfile(output_dir, clip_dir), 'dir')
%         mkdir(fullfile(output_dir, "result_clips"));
%     end
%     
%     % Keep relevant clipped variables from start of echogram for comparions
%     % of overlapping results
%     clip_start = struct('Easting', radar_tmp.Easting(1:2*clip), ...
%         'Northing', radar_tmp.Northing(1:2*clip), ...
%         'groups', radar_tmp.groups(:,1:2*clip),...
%         'likelihood', radar_tmp.likelihood(:,1:2*clip), ...
%         'depth',radar_tmp.depth, 'ages',radar_tmp.ages(:,1:2*clip,:));
%     clip_start.SMB_yr = radar_tmp.SMB_yr(1:2*clip);
%     clip_start.SMB = radar_tmp.SMB(1:2*clip);
%     
%     % Generate file names/paths for starting clip output and save
%     fn_start = sprintf('%s%d%s','clip_start',i, '.mat');
%     start_path = fullfile(output_dir, clip_dir, fn_start);
%     [save_success2] = parsave(clip_start, start_path);
%     
%     % Keep relevant clipped variables from end of echogram for comparions
%     % of overlapping results
%     clip_end = struct('Easting', radar_tmp.Easting(end-2*clip+1:end), ...
%         'Northing', radar_tmp.Northing(end-2*clip+1:end), 'groups', ...
%         radar_tmp.groups(:,end-2*clip+1:end), 'likelihood', ...
%         radar_tmp.likelihood(:,end-2*clip+1:end), 'depth', ...
%         radar_tmp.depth, 'ages', radar_tmp.ages(:,end-2*clip+1:end,:));
%     clip_end.SMB_yr = radar_tmp.SMB_yr(end-2*clip+1:end);
%     clip_end.SMB = radar_tmp.SMB(end-2*clip+1:end);
%     
%     % Generate file names/paths for end clip output and save
%     fn_end = sprintf('%s%d%s','clip_end',i, '.mat');
%     end_path = fullfile(output_dir, clip_dir, fn_end);
%     [save_success3] = parsave(clip_end, end_path);
    
end

% delete(tmp_path)
