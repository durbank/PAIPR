% This file concatenates locations from all the GPR files in this folder
% in order to determine their real-world locations

% % Directory to radar files of interest
% radar_dir = ['SEAT_Traverses' filesep 'core_grids' filesep ...
%      'grid_SEAT10_1' filesep];

function [locations] = concat_loc(file_dir, wild)

% List all files matching 'wild' within radar directory
files = dir(strcat(file_dir, wild));

lat0 = [];
lon0 = [];
for file = files'
    load(strcat(file_dir, file.name), 'lat', 'lon')
    lat0 = [lat0 lat];
    lon0 = [lon0 lon];
end

lat = lat0;
lon = lon0;
idx = logical(lat <= -60);
clear file files lat0 lon0

% Convert from lat/lon to polar coordinates (requires the Antarctic Mapping
% Toolbox (AMT)
[Easting, Northing] = ll2ps(lat(idx), lon(idx));

locations = struct('lat', lat, 'lon', lon, 'lat_idx', idx, 'Easting',...
    Easting, 'Northing', Northing);

