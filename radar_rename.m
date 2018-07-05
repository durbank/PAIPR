% Small script to rename radar matlab files (this is sometimes neccessary
% to make naming conventions consistent across all files in order to
% properly import and process data in the correct order). This script
% copies files to a new directory, leaving the original data (with original
% naming scheme) intact


% Directories to data of interest based on computer
PC_true = ispc;
switch PC_true
    case true
        data_path = 'F:/Research/Antarctica/Data/';
        addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Addons/';
end

% Paths to the source and destination data directories
radar_dir = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband/SEAT10_4toSEAT10_6');
parts = strsplit(radar_dir, '/');
parent_path = strjoin(parts(1:end-1), '/');
output_dir = fullfile(parent_path, strcat(parts{end}, '_new/'));

mkdir(output_dir);


% Get list of radar files in directory
wild = '*.mat';
files = dir(fullfile(radar_dir, wild));

for i = 1:length(files)
    if contains(files(i).name, 'Ku_band') == true
        fn = strcat('layers_ku_band_', files(i).name(16:end));
    else
        fn = files(i).name;
    end
    output_path = fullfile(output_dir, fn);
    copyfile(fullfile(files(i).folder, files(i).name), output_path);
end