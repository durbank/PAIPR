% Function to take a long-formatted .csv file of density profiles, nest on
% unique locations, and convert lat/lon to Easting/Northing

function [rho_table] = import_rho(csv_file)

T_long = readtable(csv_file);
locs = unique(T_long(:,{'Lat', 'Lon'}), 'rows');


var_Data = cell(height(locs), 1);
for i=1:length(var_Data)
    i_idx = T_long.Lat == locs.Lat(i) & T_long.Lon == locs.Lon(i);
    var_Data{i} = T_long(i_idx, {'Depth', 'pred_mean', 'pred_sd'});
end

% Determine whether data is in Greenland or Antarctica
if mean(locs{:,'Lat'}) > 0
    % Generate projection structure for EPSG:3413 (NSIDC Sea Ice Polar
    % Stereographic North)
    proj = defaultm('ups');
    proj.geoid = wgs84Ellipsoid('meters');
    proj.maplatlimit = [84, 90];
    proj.maplonlimit = [-180, 180];
    proj.origin = [90,0,0];
    proj.flatlimit = [-Inf,6];
    proj.flonlimit = [-180,180];
    
    % Convert lat/lon coordinates to Easting/Northing
    [rho_E, rho_N] = projfwd(proj,locs{:,'Lat'},locs{:,'Lon'});
    
else
    % Convert lat/lon coordinates to Easting/Northing
    [rho_E, rho_N] = ll2ps(locs{:,'Lat'}, locs{:,'Lon'});
end


rho_table = table(rho_E, rho_N, var_Data, 'VariableNames', ...
    {'Easting', 'Northing', 'Data'});

end