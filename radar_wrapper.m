% Wrapper function to import and process an arbitrary number of radargram
% files within a single directory

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'laptop';
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
radar_dir = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband/', ...
    'SEAT10_4toSEAT10_6');
% radar_dir = fullfile(data_path, 'IceBridge/SEAT10_4to10_6/2011_SNO');

%%

radar_ALL = radar_format(radar_dir);
overlap = 10000;
horz_res = 25;
output_dir = 'SMB_results';

for i = 1:length(radar_ALL)
    
    [radar_tmp] = radar_RT(radar_ALL(i).segment, cores, Ndraw);
    [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
    
    clip = round(0.5*overlap/horz_res);
    
    fld_nm = fieldnames(radar_tmp);
    fld_want = {'collect_date', 'Easting', 'Northing', 'dist', 'depth', ...
        'data_smooth', 'peaks', 'groups', 'ages', 'SMB_yr', 'SMB'};
    
    radar = struct('collect_date', radar_tmp.collect_date, ...
        'Easting', radar_tmp.Easting(clip:end-clip),...
        'Northing', radar_tmp.Northing(clip:end-clip), ...
        'dist', radar_tmp.dist(clip:end-clip), 'depth', radar_tmp.depth, ...
        'data_smooth', radar_tmp.data_smooth(:,clip:end-clip),...
        'peaks', radar_tmp.peaks(:,clip:end-clip), ...
        'groups', radar_tmp.groups(:,clip:end-clip),...
        'likelihood', radar_tmp.likelihood(:,clip:end-clip), ...
        'ages', radar_tmp.ages(:,clip:end-clip,:));
    radar.dist = radar.dist - radar.dist(1);
    if isfield(radar_tmp, 'elev')
        radar.elev = radar_tmp.elev(clip:end-clip);
    end
    radar.SMB_yr =  radar_tmp.SMB_yr(clip:end-clip);
    radar.SMB = radar_tmp.SMB(clip:end-clip);
    
    filename = sprintf('%s%d%s','radar_out',i, '.mat');
    output = fullfile(radar_dir, output_dir, filename);
    save(output, '-struct', 'radar', '-v7.3')
    
end
