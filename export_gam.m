% Tmp script to take processed PAIPR data, fit gamma distributions to the
% SMB time series, and reformat as a long form time series data table for
% easy import into R

[file, path] = uigetfile('*.mat', ...
    'Select processed PAIPR .mat file to reformat');
data = load(fullfile(path, file));

% output_dir = uigetdir('', 'Select directory to save converted tables');



% Add Antarctic Mapping Toolbox (AMT) to path
addon_path = '/home/durbank/MATLAB/Add-Ons/';
addon_folder = fullfile(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

% Convert Northing, Easting to lat/lon
[data.lat, data.lon] = ps2ll(data.Easting, data.Northing);



table_cell = cellfun(@(a,b,c,d,e) accum_distGamma(a, b, c, d, e), ...
    num2cell(data.lat), num2cell(data.lon), num2cell(data.elev), ...
    data.SMB, data.SMB_yr, 'UniformOutput', false);
long_table = vertcat(table_cell{:});


output_dir = uigetdir('', 'Select directory to save converted tables');
f_name = 'mat_table.csv';
writetable(long_table, fullfile(output_dir, f_name));





