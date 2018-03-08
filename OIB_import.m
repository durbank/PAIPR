% Function to import and format Operation IceBridge data for use with
% accum-radar functions

function [mdata] = OIB_import(file)

% Load data as a .mat structure
mdata = load_L1B(file);

%% Elevation correction
param = [];
param.update_surf = true;
param.filter_surf = true;
param.er_ice = 3.15;
param.depth = '[min(Surface_Elev)-20 max(Surface_Elev)+2]';
[mdata_c,depth_good_idxs] = elevation_compensation(mdata,param);

data_depth = size(mdata.Data, 1)-max(mdata_c.Surface_Bin);
data_out = zeros(data_depth, size(mdata.Data, 2));
TWTT = zeros(data_depth, size(mdata.Data, 2));
for i = 1:size(mdata.Data, 2)
    data_out(:,i) = mdata.Data(mdata_c.Surface_Bin(i):mdata_c.Surface_Bin(i)+data_depth-1,i);
    TWTT(:,i) = mdata.Time(mdata_c.Surface_Bin(i):mdata_c.Surface_Bin(i)+data_depth-1);
end

data_depth_min = min(sum(data_out>0));
data_out = data_out(1:data_depth_min,:);
TWTT = TWTT(1:data_depth_min,:);
TWTT = TWTT(:,1) - TWTT(1,1);

%% Cut off data below ~35 m depth

% Calculate velocity (assuming constant relative permittivity)
c = 2.9979E8;
u = c/sqrt(param.er_ice);
idx_end = ceil((35/u)/(0.5*TWTT(2)));

% Truncate data at this cut off point
data_out = data_out(1:idx_end,:);
TWTT = TWTT(1:idx_end);

% Convert struct variable names to match those for 'radar_clean.m'
lat = mdata.Latitude';
lon = mdata.Longitude';
elev = mdata.Elevation';
time_gps = mdata.GPS_time';
time_trace = mode(diff(TWTT));
data_out = 10*log10(data_out);

%% Modified code from 'radar_clean.m' (eventually want to combine this into single funtion)

% Remove data without valid location values (lat/lon)
loc_idx = logical(lat);
lat(~loc_idx) = [];
lon(~loc_idx) = [];
time_gps(~loc_idx) = [];
time_trace(~loc_idx) = [];
data_out(:,~loc_idx) = [];

% Calucate distance along traverse (in meters)
d = pathdist(lat, lon);

% Remove data from repeat traces at the same location
dist_idx = logical(diff(d));
d(~dist_idx) = [];
lat(~dist_idx) = [];
lon(~dist_idx) = [];
time_gps(~dist_idx) = [];
time_trace(~dist_idx) = [];
data_out(:,~dist_idx) = [];

% Identify data columns consisting solely of NaN cells, and remove these
% columns from all radar data arrays
nan_idx = all(isnan(data_out));
d(nan_idx) = [];
lat(nan_idx) = [];
lon(nan_idx) = [];
time_gps(nan_idx) = [];
time_trace(nan_idx) = [];
data_out(:,nan_idx) = [];

% Fill in missing radar data (missing/nan values) column-wise using a
% piecewise shape-preserving spline
data_out = fillmissing(data_out, 'pchip');

% Convert lat/lon coordinates to polar stereo coordinates (uses AMT)
[Easting, Northing] = ll2ps(lat, lon);

% Define time (decimal calendar year) at the surface (the time of data 
% collection with zero time for data defined as 1970 Jan 01 00:00:00)
[~, name, ~] = fileparts(file);
date_char = char(extractBetween(name,'_','_'));
collect_date = datetime(str2double(date_char(1:4)), ...
    str2double(date_char(5:6)), str2double(date_char(7:8)));
% collect_date = datetime(1970,1,1,0,0,0) + seconds(mean(mdata.time_gps));
collect_date = year(collect_date) + day(collect_date, 'dayofyear')/365;

% Create new mdata structure from requisite variables
mdata = struct('lat', lat, 'lon', lon, 'Northing', Northing, 'Easting', ...
    Easting,  'dist', d, 'elev', elev, 'data_out', data_out, 'TWTT', ...
    TWTT, 'time_trace', time_trace, 'collect_date', collect_date);

end
