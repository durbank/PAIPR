% Function to import and format Operation IceBridge data for use with
% accum-radar functions

% Varargin defines the starting and ending trace indices to include in the
% data import. These should numeric scalars, with values between 1 and
% the total length of the radargram file to import. If both starting and
% ending indices are prescribed, varargin{2} should be greater than
% varargin{1}

function [mdata] = OIB_convert(file)

% Load data as a .mat structure
mdata = load_L1B(file);

%% Elevation correction
param = [];
param.update_surf = true;
param.filter_surf = true;
param.er_ice = 3.15;
param.depth = '[min(Surface_Elev)-20 max(Surface_Elev)+2]';
[mdata_c, ~] = elevation_compensation(mdata,param);

% Interpolate missing values in mdata_c-Surface_bin
mdata_c.Surface_Bin = round(fillmissing(mdata_c.Surface_Bin, 'spline'));

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

%% Cut off data below ~30 m depth

% Calculate velocity (assuming constant relative permittivity)
c = 2.9979E8;
u = c/sqrt(2.25);
depth_cut = 30;     % Approximate cut-off depth (data below this depth unreliable)
idx_end = min([ceil((depth_cut/u)/(0.5*TWTT(2))) size(data_out,1)]);

% Truncate data at this cut off point
data_out = data_out(1:idx_end,:);
TWTT = TWTT(1:idx_end);

% Get collection time for each trace in echogram
ncid = netcdf.open(file);
time_id = netcdf.inqVarID(ncid, 'time');
start_char = netcdf.getAtt(ncid, time_id, 'units');
char_split = split(start_char, 'seconds since ');
t_start = datetime(strcat(char_split{2}), ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss');
collect_time = t_start + seconds(mdata.GPS_time);
netcdf.close(ncid);

% % Define time (decimal calendar year) at the surface (the time of data 
% % collection with zero time for data defined as 1970 Jan 01 00:00:00)
% [~, name, ~] = fileparts(file);
% date_char = char(extractBetween(name,'_','_'));
% collect_date = datetime(str2double(date_char(1:4)), ...
%     str2double(date_char(5:6)), str2double(date_char(7:8)));
% collect_date = year(collect_date) + day(collect_date, 'dayofyear')/365;

% Create new mdata structure from existing variables such that the naming
% conventions match those from SEAT radar fieldnames
mdata = struct('lat', mdata.Latitude', 'lon', mdata.Longitude', 'elev', ...
    mdata.Elevation', 'data_out', 10*log10(data_out), 'collect_time', ...
    collect_time', 'time_trace', repmat(0.5*mode(diff(TWTT)), 1, ...
    length(mdata.Latitude)));
% mdata = struct('lat', mdata.Latitude', 'lon', mdata.Longitude', 'elev', ...
%     mdata.Elevation', 'data_out', 10*log10(data_out), 'time_gps', ...
%     mdata.GPS_time', 'time_trace', repmat(0.5*mode(diff(TWTT)), 1, ...
%     length(mdata.Latitude)), 'collect_time', collect_time');

end
