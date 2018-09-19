% Wrapper function to import and process an arbitrary number of radargram
% files within a single directory

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'work';
        %         computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/WARP backup/Research/Antarctica/Data/';
%                 data_path = 'E:/Research/Antarctica/Data/';
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

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Define number of Monte Carlo simulations to perform
Ndraw = 100;

% Define path to the directory containing radar data (relative to the
% 'data' directory path)
% radar_dir = fullfile(data_path, 'radar/SEAT_Traverses/core-site_tests/', ...
%     'SEAT10_5');
radar_dir = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband/', ...
    'SEAT10_4toSEAT10_6');

% List all files matching 'wild' within radar directory
wild = '*.mat';
files = dir(fullfile(radar_dir, wild));

% If directory has no '.mat' files, additionally check for '.nc' files
if isempty(files)
    wild = '*.nc';
    files = dir(fullfile(radar_dir, wild));
end

%%

% Preallocate arrays for continuous lat/lon positions for all data in
% directory
lat = [];
lon = [];

% Preallocate array for the radargram length (in data bins) for each 
% data file in directory
file_length = zeros(1, length(files));

% Iteratively load lat/lon positions for each data file, and add
% corresponding values to arrays
for i = 1:length(files)
    lat_i = load(fullfile(files(i).folder, files(i).name), 'lat');
    lon_i = load(fullfile(files(i).folder, files(i).name), 'lon');
    file_length(i) = length(lat_i.lat);
    lat = [lat lat_i.lat];
    lon = [lon lon_i.lon];
end

% Calculate the cummulative radargram length (in data bins) for the data
% files in directory
file_idx = [0 cumsum(file_length)];

% Replace data without valid location values (outside Antarctic Circle) 
% with nearest preceding valid location (these missing data will later be
% removed in processing)
invld_idx = lat>=-65 | lat<-90;
lat(invld_idx) = NaN;
lon(invld_idx) = NaN;
lat = fillmissing(lat, 'previous');
lon = fillmissing(lon, 'previous');

% Calucate distance along traverse (in meters)
d = pathdist(lat, lon);

% Find indices to break radargrams based on absence of data across an
% extended distance (greater than 500 m)
break_idx = [0 find([0 diff(d)]>500) length(lat)];

% Set the minimum length needed for radargram processing and radargram 
% overlap interval (in meters)
length_min = 30000;
overlap = 5000;

% For each continuous section of radargram (no significant breaks),
% determine additional breakpoint indices based on the desired processing
% length and degree of overlap
for i = 1:length(break_idx)-1
    
    % Determine lat/lon and distance along the ith section radar data
    lat_i = lat(break_idx(i)+1:break_idx(i+1));
    lon_i = lon(break_idx(i)+1:break_idx(i+1));
    dist_i = pathdist(lat_i, lon_i);
    
    % Initialize while loop for current radar segment
    search = true;
    j_start = break_idx(i)+1;
    
    % 
    if dist_i(end) <= length_min
        search = false;
        j_end = length(dist_i);
        [j_start j_end+break_idx(i)]
%         coords{j} = [j_start j_end];
        j = j+1;
    end
    
    while search == true
        
        dist_j = dist_i - dist_i(j_start-break_idx(i));
        j_end = find(dist_j >= length_min, 1);
        
        if dist_j(end)-dist_j(j_end) <= overlap
            
            search = false;
            j_end = length(dist_j);
        end
        
        if isempty(j_end)
            
            search = false;
            j_start = find(dist_j>=dist_j(end)-length_min, 1) + break_idx(i);
            j_end = length(dist_i);
            
        end
        
        [j_start j_end+break_idx(i)]
        j_start = find(dist_j-dist_j(j_end)+overlap>=0, 1) + break_idx(i);
        j = j+1;
    end
end


