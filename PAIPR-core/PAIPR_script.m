% Wrapper function to import and process an arbitrary number of radargram
% files within a single directory

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

% Define path to the directory containing radar data (relative to the
% 'data' directory path)
radar_dir = uigetdir(data_path, ...
    'Select folder containing .mat radar data files');
% radar_dir = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband/', ...
%     'SEAT10_4toSEAT10_6');
% radar_dir = fullfile(data_path, 'IceBridge/SEAT10_4to10_6/2011_SNO');

% Define path to output directory of where to store results
output_dir = uigetdir(data_path, ...
    "Select folder to store output .mat results files");

%%

% Combine all input echograms and decompose into echogram subdomains of
% desired length
radar_ALL = echo_format(radar_dir);

% Define overlap distance and the final horizontal resolution of the output
% data
overlap = 10000;
horz_res = 25;

% Determine if all data break segments are sufficiently long to process
keep_idx = false(length(radar_ALL), 1);
for i = 1:length(radar_ALL)
   
    if radar_ALL(i).segment.dist(end) >= 1.5*overlap
        keep_idx(i) = true;
    end
end

% Only keep data segments of sufficient length
radar_ALL = radar_ALL(keep_idx);

% Parellel for loop to process all decomposed echograms
parfor i = 1:length(radar_ALL)
    
    % Calculate radar age-depth profile distributions (includes processing
    % steps for depth, signal-noise, radon transforms, layer tracing,
    % likelihood assignments, and age calculations)
    r = -4.3491e-4;
    k = 4.600;
    [radar_tmp] = calc_age(radar_ALL(i).segment, cores, r, k, Ndraw);
    
    % Calculate radar annual SMB
    [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
    
    
    %%
    
%     % May be used in future versions of code
%     fld_nm = fieldnames(radar_tmp);
%     fld_want = {'collect_date', 'Easting', 'Northing', 'dist', 'depth', ...
%         'rho_coeff', 'rho_var', 'data_smooth', 'peaks', 'groups', 'ages', ...
%         'SMB_yr', 'SMB'};
    
    % Clip radar data structure variables based on the desired radargram
    % overlap, and combine desired clipped variables into new data
    % structure
    clip = round(0.5*overlap/horz_res);
    radar = struct('collect_date', radar_tmp.collect_date, ...
        'Easting', radar_tmp.Easting(clip:end-clip),...
        'Northing', radar_tmp.Northing(clip:end-clip), ...
        'dist', radar_tmp.dist(clip:end-clip), 'depth', radar_tmp.depth, ...
        'rho_coeff', radar_tmp.rho_coeff(:,clip:end-clip), ...
        'rho_var', radar_tmp.rho_var(:,clip:end-clip),...
        'data_smooth', radar_tmp.data_smooth(:,clip:end-clip),...
        'peaks', radar_tmp.peaks(:,clip:end-clip), ...
        'groups', radar_tmp.groups(:,clip:end-clip),...
        'likelihood', radar_tmp.likelihood(:,clip:end-clip), ...
        'ages', radar_tmp.ages(:,clip:end-clip,:));
    
    % Redefine radargram distances based on clipped data
    radar.dist = radar.dist - radar.dist(1);
    
    % Check for existing elevation data, and add to output struct if
    % available
    if isfield(radar_tmp, 'elev')
        radar.elev = radar_tmp.elev(clip:end-clip);
    end
    
    % Clip SMB cell array data and add to output struct
    radar.SMB_yr =  radar_tmp.SMB_yr(clip:end-clip);
    radar.SMB = radar_tmp.SMB(clip:end-clip);
    
    % Generate file names and paths under which to save data
    filename = sprintf('%s%d%s','radar_out',i, '.mat');
    output = fullfile(output_dir, filename);
    
    % Save output structures to disk
    [save_success1] = parsave(radar, output)
    
    %% Output overlapping data for comparisons of neighboring echograms
    
    % Check for existence of directory for clipped results in output_dir,
    % and create it if it does not exist
    clip_dir = "result_clips";
    if ~exist(fullfile(output_dir, clip_dir), 'dir')
        mkdir(fullfile(output_dir, "result_clips"));
    end
    
    % Keep relevant clipped variables from start of echogram for comparions 
    % of overlapping results
    clip_start = struct('Easting', radar_tmp.Easting(1:2*clip), ...
        'Northing', radar_tmp.Northing(1:2*clip), ...
        'groups', radar_tmp.groups(:,1:2*clip),...
        'likelihood', radar_tmp.likelihood(:,1:2*clip), ...
        'depth',radar_tmp.depth, 'ages',radar_tmp.ages(:,1:2*clip,:));
    clip_start.SMB_yr = radar_tmp.SMB_yr(1:2*clip);
    clip_start.SMB = radar_tmp.SMB(1:2*clip);
    
    % Generate file names/paths for starting clip output and save
    fn_start = sprintf('%s%d%s','clip_start',i, '.mat');
    start_path = fullfile(output_dir, clip_dir, fn_start);
    [save_success2] = parsave(clip_start, start_path)
    
    % Keep relevant clipped variables from end of echogram for comparions 
    % of overlapping results
    clip_end = struct('Easting', radar_tmp.Easting(end-2*clip+1:end), ...
        'Northing', radar_tmp.Northing(end-2*clip+1:end), ...
        'groups', radar_tmp.groups(:,end-2*clip+1:end),...
        'likelihood', radar_tmp.likelihood(:,end-2*clip+1:end), ...
        'depth',radar_tmp.depth, 'ages',radar_tmp.ages(:,end-2*clip+1:end,:));
    clip_end.SMB_yr = radar_tmp.SMB_yr(end-2*clip+1:end);
    clip_end.SMB = radar_tmp.SMB(end-2*clip+1:end);
    
    % Generate file names/paths for end clip output and save
    fn_end = sprintf('%s%d%s','clip_end',i, '.mat');
    end_path = fullfile(output_dir, clip_dir, fn_end);
    [save_success3] = parsave(clip_end, end_path)
    
end
