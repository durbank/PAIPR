% Function to import estimated depth-density curves and density variance at
% trace locations from previously-generated results using the stats model
% with monotonicity constraints (in the future this function should run the
% requisite R script to extract the density curves at the specified
% echogram locations using the previously-optimized stats model)

function [trace_rho] = load_rho(rho_file, Easting, Northing)

% Read in density data from external file
rho_mod = readtable(rho_file);

% Convert Easting/Northing to lat/lon;
[E, N] = ll2ps(rho_mod{:,'Lat'}, rho_mod{:,'Lon'});
T_locs = table(E,N,'VariableNames', {'Easting', 'Northing'});
rho_convert = [T_locs rho_mod(:,3:end)];

% Find unique locations within loaded density file
rho_loc = unique(rho_convert(:,{'Easting', 'Northing'}));

% Find nearest depth-density profile to each trace in echogram
near_idx = zeros(1, length(Easting));
for i = 1:length(Easting)
    
    dist = (Easting(i) - rho_loc.Easting).^2 + ...
        (Northing(i) - rho_loc.Northing).^2;
    [~, near_idx(i)] = min(dist);
        
end
trace_rho = rho_loc(near_idx,:);


% Preallocate cell array for depth-density data at each location
trace_rho.Data = cell(height(trace_rho),1);

% Add depth-density data to table for each location
for i = 1:height(trace_rho)
    
    data_idx = rho_convert.Easting == trace_rho.Easting(i) & ...
        rho_convert.Northing == trace_rho.Northing(i);
    trace_rho.Data{i} = rho_convert(data_idx,3:end);
end






end