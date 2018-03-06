% Script to import and format Operation IceBridge radar data for use with
% accum-radar scripts

PC_true = ispc;
switch PC_true
    case true
        data_path = 'D:/Research/Antarctica/WAIS Variability/';
        addon_path = 'D:/Research/Antarctica/WAIS Variability/Addons/';
    case false
        data_path = '/Volumes/WARP/Research/Antarctica/WAIS Variability/';
        addon_path = '/Users/Durbank/Documents/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

addpath cresis-L1B-matlab-readers/

% Files IRSNO1B_20111109_02_211 through ...02_227 roughly follows the
% SEAT10-1 to SEAT10-6 transect. The same holds true for the following:
%   IRSNO1B_20161109_02_320 through IRSNO1B_20161109_02_336
%   Ku-band: IRKUB1B_20161109_02_320 through IRKUB1B_20161109_02_336

% Files IRSNO1B_20111109_02_242 through ...02_272 exactly follows the 
% SEAT10-4 - SEAT10-6 transect. The same holds true for the following:
%   IRSNO1B_20161109_02_350 through IRSNO1B_20161109_02_381
%   Ku-band: IRKUB1B_20161109_02_350 through IRKUB1B_20161109_02_381

file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Kuband/2016/IRKUB1B_20161109_02_381.nc';

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
TWTT = TWTT - TWTT(1,:);


% Convert struct variable names to match those for 'radar_clean.m'
mdata.lat = mdata.Latitude';
mdata.lon = mdata.Longitude';
mdata.time_gps = mdata.GPS_time';
mdata.time_trace = mode(diff(TWTT));
mdata.data_out = data_out;

%% Modified code from 'radar_clean.m' (eventually want to combine this into single funtion)

% Remove data without valid location values (lat/lon)
loc_idx = logical(mdata.lat);
mdata.lat(~loc_idx) = [];
mdata.lon(~loc_idx) = [];
mdata.time_gps(~loc_idx) = [];
mdata.time_trace(~loc_idx) = [];
mdata.data_out(:,~loc_idx) = [];

% Calucate distance along traverse (in meters)
mdata.dist = pathdist(mdata.lat, mdata.lon);

% Remove data from repeat traces at the same location
dist_idx = logical(diff(mdata.dist));
mdata.dist(~dist_idx) = [];
mdata.lat(~dist_idx) = [];
mdata.lon(~dist_idx) = [];
mdata.time_gps(~dist_idx) = [];
mdata.time_trace(~dist_idx) = [];
mdata.data_out(:,~dist_idx) = [];

% Identify data columns consisting solely of NaN cells, and remove these
% columns from all radar data arrays
nan_idx = all(isnan(mdata.data_out));
mdata.dist(nan_idx) = [];
mdata.lat(nan_idx) = [];
mdata.lon(nan_idx) = [];
mdata.time_gps(nan_idx) = [];
mdata.time_trace(nan_idx) = [];
mdata.data_out(:,nan_idx) = [];

% Fill in missing radar data (missing/nan values) column-wise using a
% piecewise shape-preserving spline
mdata.data_out = fillmissing(mdata.data_out, 'pchip');

% Convert lat/lon coordinates to polar stereo coordinates (uses AMT)
[mdata.Easting, mdata.Northing] = ll2ps(mdata.lat, mdata.lon);

% Define time of the surface
mdata.collect_date = mean(mdata.time_gps);

% Remove uneeded fields from the mdata structure
mdata = rmfield(mdata, {'Latitude', 'Longitude', 'GPS_time', 'time_gps'});