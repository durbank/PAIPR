% Script to import and clean GPR profiles, and prepare it for further
% analysis

% Varargin defines the starting and ending trace indices to include in the
% data import. These should be numeric scalars, with values between 1 and
% the total length of the radargram file to import. If both starting and
% ending indices are prescribed, varargin{2} should be greater than
% varargin{1}

function [mdata] = import_radar(file, varargin)

% Load the radar data file
if isstruct(file)
    mdata = file;
else
    mdata = load(file);
end

% Determine trace starting and ending indices from varargin variable
if isempty(varargin)
    start_idx = 1;
    end_idx = length(mdata.lat);
elseif length(varargin) == 1
    start_idx = varargin{1};
    end_idx = length(mdata.lat);
elseif isempty(varargin{1})
    start_idx = 1;
    end_idx = varargin{2};
else
    start_idx = varargin{1};
    end_idx = varargin{2};
end

% Calucate distance along traverse (in meters)
distances = pathdist(mdata.lat, mdata.lon);
        
bounds_idx = 1:length(mdata.lat) < start_idx | 1:length(mdata.lat) > end_idx;
loc_idx = mdata.lat >= -65;
dist_idx = ~logical([1 diff(distances)]);
nan_idx = all(isnan(mdata.data_out));

rm_idx = logical(sum([bounds_idx; loc_idx; dist_idx; nan_idx]));

mdata.lat(rm_idx) = [];
mdata.lon(rm_idx) = [];
mdata.time_gps(rm_idx) = [];
mdata.time_trace(rm_idx) = [];
mdata.data_out(:,rm_idx) = [];

if ~isstruct(file)
    
    % Approximate the surface in the radar trace, and remove data above
    % the surface (this could be improved in the future. OIB may have
    % some code that could be helpful)
    surf_row = 203;
    air_idx = 1:size(mdata.data_out,1) < surf_row;
    mdata.data_out(air_idx,:) = [];
    
    if isfield(mdata, 'arr_layers')
        mdata.arr_layers(:,rm_idx) = [];
        mdata.arr_layers(air_idx,:) = [];
    end
    
end

if isfield(mdata, 'elev')
    mdata.elev(rm_idx) = [];
end



% Fill in missing radar data (missing/nan values) column-wise using a
% piecewise shape-preserving spline
mdata.data_out = fillmissing(mdata.data_out, 'pchip');

% Convert lat/lon coordinates to polar stereo coordinates (uses AMT)
[mdata.Easting, mdata.Northing] = ll2ps(mdata.lat, mdata.lon);

% Calucate path distance along traverse (in meters)
try
    mdata.dist = pathdist(mdata.lat, mdata.lon);
catch
    mdata.dist = [];
    disp('Warning: error with calculating pathdist (Line 81 of import_radar.m)')
end

% Check to see if collection data is already present. If it is not, create
% collect_date field
if ~isfield(mdata, 'collect_date')
    % Define the age of the top of the radar (the data collection date)
    % (this will require modification when incorporating data beyond SEAT
    % traverses)
    if contains(file, '2011') == true
        DateString = '15.12.2011';
        vector = datevec(DateString, 'dd.mm.yyyy');
        day_of_year = datenum(vector) - datenum(vector(1), 1, 1);
        mdata.collect_date = vector(1) + day_of_year/365.25;
    elseif contains(file, '2010') == true
        DateString = '15.12.2010';
        vector = datevec(DateString, 'dd.mm.yyyy');
        day_of_year = datenum(vector) - datenum(vector(1), 1, 1);
        mdata.collect_date = vector(1) + day_of_year/365.25;
    else
        disp('Check collection date')
    end
end

end