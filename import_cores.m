% Function to import core data from .xlsx file (specific format required)
% and arrange it in a structure for future MATLAB analysis

% core_file = '/Volumes/WARP DRIVE/Research/Antarctica/WAIS Variability/SEAT_cores/core_data_interp.xlsx';
function [cores] = import_cores(core_file)

% % Add Antarctic Mapping Toolbox (AMT) to path
% addon_path = ['Addons' filesep 'AntarcticMappingTools_v5.03' filesep];
% addpath(genpath(addon_path))

% Import basic data about firn cores in .xlsx file (name, location,
% elevation, etc) and places inside a data structure
[data, name] = xlsread(core_file, 'Site Information');
name = name(1,2:end);
cores = struct('name', {name}, 'lat', data(1,:), 'lon', data(2,:), 'elev', ...
    data(3,:));

% Convert core lat/lon locations to polar sterographic projection (meters)
[cores.Easting, cores.Northing] = ll2ps(cores.lat, cores.lon);

% Import additional core information (depth, age, density, isotopes, etc), 
% fill in NaN/missing values in density and isotope data, and place data in 
% the 'cores' data structure
for i = 1:numel(cores.name)
    [data, ~] = xlsread(core_file, cores.name{i});
    cores.(cores.name{i}) = struct('name', cores.name(i), 'Easting', ...
        cores.Easting(i), 'Northing', cores.Northing(i), 'depth', ...
        data(:,1), 'age', data(:,2), 'rho', fillmissing(data(:,3), 'pchip'),...
        'dD', fillmissing(data(:,4), 'pchip'), ...
        'd18O', fillmissing(data(:,5), 'pchip'));
end

end