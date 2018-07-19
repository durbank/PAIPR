% Small script to contcatonate all radar files within a given directory to
% a single radar file

% THIS IS THE SCRIPT USED TO CONCATONATE 2010 SEAT RADAR ALONG THE
% SEAT10-4 TO SEAT10-6 TRAVERSE

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
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

data_name = 'SEAT10_4toSEAT10_6';
radar_dir = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband', data_name);
file_name = strcat('layers_ku_band_',data_name, '.mat');
output_dir = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband');

% Get list of radar files in directory
wild = '*.mat';
files = dir(fullfile(radar_dir, wild));

% Initialize for loop with first file in directory
data_struct = load(fullfile(radar_dir, files(1).name));
fld_names = fieldnames(data_struct);
data_cells = struct2cell(data_struct);
data_cells = data_cells(1:5);

for i = 2:length(files)
    data_i = struct2cell(load(fullfile(radar_dir, files(i).name)));
    data_cells = cellfun(@horzcat, data_cells, data_i(1:5), 'UniformOutput', false);
end

output_struct = cell2struct(data_cells, fld_names(1:5), 1);

% lat = output_struct.lat;
% lon = output_struct.lon;
% time_gps = output_struct.time_gps;
% time_trace = output_struct.time_trace;
% data_out = output_struct.data_out;
% 
% Define the age of the top of the radar (the date the radar was collected)
% (this will require modification when incorporating data beyond SEAT
% traverses)
if contains(files(i).name, '2011') == true
    DateString = '01.01.2012';
    vector = datevec(DateString, 'dd.mm.yyyy');
    day_of_year = datenum(vector) - datenum(vector(1), 1, 1);
    collect_date = vector(1) + day_of_year/365.25;
elseif contains(files(i).name, '2010') == true
    DateString = '01.01.2011';
    vector = datevec(DateString, 'dd.mm.yyyy');
    day_of_year = datenum(vector) - datenum(vector(1), 1, 1);
    collect_date = vector(1) + day_of_year/365.25;
else
    disp('Check collection date')
end


output_struct.collect_date = collect_date;
output_path = fullfile(output_dir, file_name);
save(output_path, '-struct', 'output_struct', '-v7.3')

% save(output_path, 'lat', 'lon', 'time_gps', ...
%     'time_trace', 'data_out', 'collect_date', '-v7.3')

