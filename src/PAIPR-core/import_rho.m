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

[rho_E, rho_N] = ll2ps(locs{:,'Lat'}, locs{:,'Lon'});
rho_table = table(rho_E, rho_N, var_Data, 'VariableNames', ...
    {'Easting', 'Northing', 'Data'});

end