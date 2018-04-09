% Small script to contcatonate all radar files within a given directory to
% a single radar file

% % User defines radar directory path
% radar_dir = input('Path to radar directory:', 's');
% if radar_dir(end) ~= filesep
%     radar_dir = strcat(radar_dir, filesep);
% end
% 
% % User defines radar output location
% output_dir = input('Location to save output to:', 's');
% if output_dir(end) ~= filesep
%     output = strcat(output_dir, filesep);
% end

radar_dir = 'E:\Research\Antarctica\WAIS Variability\SEAT_Traverses\SEAT2010Kuband\ProcessedSEAT2010\transectSEAT10_5_6\';
output_dir = 'E:\Research\Antarctica\Data\OUTPUT\';
output_name = 'layers_ku_band_transectSEAT10_5_6.mat';

% Get list of radar files in directory
wild = 'layers*';
files = dir(strcat(radar_dir, wild));

% Initialize for loop with first file in directory
data_struct = load(strcat(radar_dir, files(1).name));
fld_names = fieldnames(data_struct);
data_cells = struct2cell(data_struct);
data_cells = data_cells(1:5);

for i = 2:length(files)
    data_i = struct2cell(load(strcat(radar_dir, files(i).name)));
    data_cells = cellfun(@horzcat, data_cells, data_i(1:5), 'UniformOutput', 0);
end

output_struct = cell2struct(data_cells, fld_names(1:5), 1);

lat = output_struct.lat;
lon = output_struct.lon;
time_gps = output_struct.time_gps;
time_trace = output_struct.time_trace;
data_out = output_struct.data_out;

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

save(strcat(output_dir, output_name), 'lat', 'lon', 'time_gps', ...
    'time_trace', 'data_out', 'collect_date', '-v7.3')

