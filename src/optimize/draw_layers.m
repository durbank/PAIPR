% Re-usable script to manually draw layers from radar echograms (using the
% 'draw_manual' function)

DATA_DIR = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'optimization/v0.4.0/flights/flight_20161109/'];
fn = 'Site_B.mat';

% File containing interim data (processed with PAIPR at least through
% 'calc_layers' function)
data_path = fullfile(DATA_DIR, 'interim_data', fn);

% Path for output of saved manully traced layers
output_path = fullfile(DATA_DIR, 'man_layers', fn);

% Load radar to workspace
radar = load(data_path);

% Run drawing routine
[man_layers] = draw_manual(radar);

%% Save manual layer output to disk for later use

% Save processed radar structure for future use
save(output_path, 'man_layers')