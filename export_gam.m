% Tmp script to take processed PAIPR data, fit gamma distributions to the
% SMB time series, and reformat as a long form time series data table for
% easy import into R

[file, path] = uigetfile('*.mat', ...
    'Select processed PAIPR .mat file to reformat');
data = load(fullfile(path, file));

% output_dir = uigetdir('', 'Select directory to save converted tables');

table_cell = cellfun(@(a,b,c,d,e) accum_distGamma(a, b, c, d, e), ...
    num2cell(data.Easting), num2cell(data.Northing), num2cell(data.elev), ...
    data.SMB, data.SMB_yr, 'UniformOutput', false);
long_table = vertcat(table_cell{:});



