% Small script to contcatonate all radar files within a given directory to
% a single radar file

% THIS IS THE SCRIPT USED TO CONCATONATE 2011 SNO OIB RADAR ALONG THE
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

% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

data_name = '2011_SNO';
radar_dir = fullfile(data_path, 'IceBridge/SEAT10_4to10_6/', data_name);
file_name = strcat(data_name, '_all.mat');
output_dir = fullfile(data_path, 'IceBridge/SEAT10_4to10_6/');

% Get list of radar files in directory
wild = '*.nc';
files = dir(fullfile(radar_dir, wild));

% Initialize for loop with first file in directory
data_struct = OIB_import(fullfile(radar_dir, files(1).name));
fld_names = fieldnames(data_struct);
data_cells = struct2cell(data_struct);

for i = 2:length(files)
    data_i = struct2cell(OIB_import(fullfile(radar_dir, files(i).name)));
    data_cells = cellfun(@horzcat, data_cells, data_i, 'UniformOutput', 0);
end

% Average values for TWTT, time_trace, and collect_date
data_cells{8} = median(data_cells{8}, 2);
data_cells{9} = median(data_cells{9});
data_cells{10} = median(data_cells{10});

% Flip data (so that the processing direction matches SEAT direction from
% SEAT10-4 to SEAT10-6)
data_flip = cellfun(@fliplr, data_cells, 'UniformOutput', false);

output_struct = cell2struct(data_flip, fld_names, 1);
output_struct.dist = pathdist(output_struct.lat, output_struct.lon);
% OIB = output_struct;

output_path = fullfile(output_dir, file_name);

% lat = output_struct.lat;
% lon = output_struct.lon;
% Northing = output_struct.Northing;
% Easting = output_struct.Easting;
% dist = output_struct.dist;
% elev = output_struct.elev;
% data_out = output_struct.data_out;
% TWTT = output_struct.TWTT;
% time_trace = output_struct.time_trace;
% collect_date = output_struct.collect_date;

% save(output_path, 'lat', 'lon', 'Northing', 'Easting',...
%     'dist', 'elev', 'data_out', 'TWTT', 'time_trace', ...
%     'collect_date', '-v7.3')

save(output_path, '-struct', 'output_struct', '-v7.3')