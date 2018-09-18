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

% List all files matching 'wild' within radar directory (either '.mat' or
% '.nc' files)
wild = '*.mat';
files = dir(fullfile(radar_dir, wild));
if isempty(files)
    wild = '*.nc';
    files = dir(fullfile(radar_dir, wild));
end

%%
lat = [];
lon = [];
file_length = zeros(1, length(files));

for i = 1:length(files)
    lat_i = load(fullfile(files(i).folder, files(i).name), 'lat');
    lon_i = load(fullfile(files(i).folder, files(i).name), 'lon');
    file_length(i) = length(lat_i.lat);
    lat = [lat lat_i.lat];
    lon = [lon lon_i.lon];
end
file_idx = [0 cumsum(file_length)];

% Replace data without valid location values (outside Arctic Circle) with
% preceding valid location
invld_idx = lat>=-65;
lat(invld_idx) = NaN;
lon(invld_idx) = NaN;
lat = fillmissing(lat, 'previous');
lon = fillmissing(lon, 'previous');

% Calucate distance along traverse (in meters)
d = pathdist(lat, lon);

% Find indices where to break radargrams based on absence of data
% across an extended difference (greater than 500 m)
break_idx = [1 find([0 diff(d)]>500) length(lat)];

% Set the minimum length needed for radargram processing
length_min = 30000;


for i = 1:length(break_idx)-1
    % Define distances along ith segment of data (definded by missing data
    % sections
    dist_i = pathdist(lat(break_idx(i):break_idx(i+1)), ...
        lon(break_idx(i):break_idx(i+1)));
    
    % Find the index where minimum length is reached
    break_i = find(dist_i>length_min, 1);
    
    
end


