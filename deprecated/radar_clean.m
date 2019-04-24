% Script to import and clean GPR profiles, and prepare it for further
% analysis
function [mdata] = radar_clean(file)

mdata = load(file);

manual_layers = isfield(mdata, 'arr_layers');

switch manual_layers
    case true
        
        % Remove data without valid location values (lat/lon)
        loc_idx = logical(mdata.lat);
        mdata.lat(~loc_idx) = [];
        mdata.lon(~loc_idx) = [];
        mdata.time_gps(~loc_idx) = [];
        mdata.time_trace(~loc_idx) = [];
        mdata.data_out(:,~loc_idx) = [];
        mdata.arr_layers(:,~loc_idx) = [];
        mdata.arr_segs(:,~loc_idx) = [];
        
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
        mdata.arr_layers(:,~dist_idx) = [];
        mdata.arr_segs(:,~dist_idx) = [];
        
        % Approximate the surface in the radar trace, and remove data above the
        % surface (this will need to be improved in the future. OIB may have some
        % code that could be helpful)
        % min_loc = 100;      % Mininum starting row value to search for surface
        % max_loc = 300;
        % [~, surf_idx] = max(abs(diff(mdata.data_out(min_loc:max_loc,:))));
        % surf_idx = surf_idx + (min_loc-1);
        % surf_row = median(surf_idx);
        % surf_row = mode(surf_idx);
        surf_row = 203;
        mdata.data_out(1:surf_row-1,:) = [];
        mdata.arr_layers(1:surf_row-1,:) = [];
        mdata.arr_segs(1:surf_row-1,:) = [];
        
        % Identify data columns consisting solely of NaN cells, and remove these
        % columns from all radar data arrays
        nan_idx = all(isnan(mdata.data_out));
        mdata.dist(nan_idx) = [];
        mdata.lat(nan_idx) = [];
        mdata.lon(nan_idx) = [];
        mdata.time_gps(nan_idx) = [];
        mdata.time_trace(nan_idx) = [];
        mdata.data_out(:,nan_idx) = [];
        mdata.arr_layers(:,nan_idx) = [];
        mdata.arr_segs(:,nan_idx) = [];
        
    case false
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
        
        % Approximate the surface in the radar trace, and remove data above the
        % surface (this will need to be improved in the future. OIB may have some
        % code that could be helpful)
        % min_loc = 100;      % Mininum starting row value to search for surface
        % max_loc = 300;
        % [~, surf_idx] = max(abs(diff(mdata.data_out(min_loc:max_loc,:))));
        % surf_idx = surf_idx + (min_loc-1);
        % surf_row = median(surf_idx);
        % surf_row = mode(surf_idx);
        surf_row = 203;
        mdata.data_out(1:surf_row-1,:) = [];
        
        % Identify data columns consisting solely of NaN cells, and remove these
        % columns from all radar data arrays
        nan_idx = all(isnan(mdata.data_out));
        mdata.dist(nan_idx) = [];
        mdata.lat(nan_idx) = [];
        mdata.lon(nan_idx) = [];
        mdata.time_gps(nan_idx) = [];
        mdata.time_trace(nan_idx) = [];
        mdata.data_out(:,nan_idx) = [];
        
end

% Fill in missing radar data (missing/nan values) column-wise using a
% piecewise shape-preserving spline
mdata.data_out = fillmissing(mdata.data_out, 'pchip');

% Convert lat/lon coordinates to polar stereo coordinates (uses AMT)
[mdata.Easting, mdata.Northing] = ll2ps(mdata.lat, mdata.lon);

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