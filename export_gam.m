% Tmp script to take processed PAIPR data, fit gamma distributions to the
% SMB time series, and reformat as a long form time series data table for
% easy import into R


% [file, path] = uigetfile('*.mat', ...
%     'Select processed PAIPR .mat file to reformat');
% data = load(fullfile(path, file));
input_dir = uigetdir('', ...
    "Select directory containing processed PAIPR .mat files to reformat");
output_dir = uigetdir('', 'Select directory to save converted tables');
files = dir(fullfile(input_dir, '*.mat'));

% Add Antarctic Mapping Toolbox (AMT) to path
addon_path = '/home/durbank/MATLAB/Add-Ons/';
addon_folder = fullfile(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

for i = 1:length(files)
    
    % Iteratively load .mat file of processed PAIPR results
    data = load(fullfile(files(i).folder, files(i).name));
    
    % Convert Northing, Easting to lat/lon
    [data.lat, data.lon] = ps2ll(data.Easting, data.Northing);
    
    % Fit gamma distribution parameters to each year of accumulation data for
    % each trace in data
    table_cell = cellfun(@(a,b,c,d,e) accum_distGamma(a, b, c, d, e), ...
        num2cell(data.lat), num2cell(data.lon), num2cell(data.elev), ...
        data.SMB, data.SMB_yr, 'UniformOutput', false);
    
    % Concatenate results into a single table in long format
    long_table = vertcat(table_cell{:});
    
    % Save the file
    [~, f_name, ~] = fileparts(fullfile(files(i).folder, files(i).name));
    writetable(long_table, fullfile(output_dir, strcat(f_name, '.csv')));
    
end





