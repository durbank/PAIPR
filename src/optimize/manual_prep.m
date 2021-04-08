% Function to process OIB Snow echograms in preparation for manually
% tracing annual layers (for PAIPR validation purposes)

function [radar] = manual_prep(radar_dir, rho_table)


% Get grouping indices of echograms
[files, start_idx, end_idx] = OIB_chunk(radar_dir);

% Define horizontal resolution for intermediate processing (in meters)
horz_res = 25;

%%

try
    echo_0 = orderfields(import_radar(...
        OIB_convert(fullfile(files(1).folder, files(1).name))));
    fields = fieldnames(echo_0);
    for i = 1:length(fields)
        echo_0.(fields{i}) = [];
    end
catch
    disp('Initial stucture assignment failed')
    disp('Reverting to manual assignment')
    
    echo_0 = struct('Easting', [], 'Northing', [], 'collect_time', [], ...
        'data_out', [], 'dist', [], 'elev', [], ...
        'lat', [], 'lon', [], 'time_trace', []);
end

for j = start_idx(1):end_idx(1)
    
    OIB_j = OIB_convert(fullfile(files(j).folder, files(j).name));
    struct_j = orderfields(import_radar(OIB_j));
    
    echo_0 = cell2struct(cellfun(@horzcat, struct2cell(echo_0), ...
        struct2cell(struct_j), 'uni', 0), fieldnames(echo_0), 1);
end

% Determine whether data are in Greenland or Antarctica
if mean(echo_0.lat) < 0
    % Calucate distance along traverse (in meters)
    echo_0.dist = pathdistps(echo_0.lat, echo_0.lon);
    
else
    % Calculate distance along traverse (in meters)
    echo_0.dist = pathdist(echo_0.lat, echo_0.lon);
    
end



[radar] = radar_stack(echo_0, horz_res);

% Load modeled depth-density data from stats model output at
% specified Easting/Northing coordinates
[rho_data] = load_rho(...
    rho_table, radar.Easting, radar.Northing);

% Convert to depth
[radar] = radar_depth(radar, rho_data);

% Calculate radar age-depth profile distributions (includes 
% processing signal-noise, radon transforms, layer tracing, 
% likelihood assignments, and age calculations)
[radar] = calc_layers(radar, 'stream');


end