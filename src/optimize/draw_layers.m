% Re-usable script to manually draw layers from radar echograms (using the
% 'draw_manual' function)

% File containing interim data (processed with PAIPR at least through
% 'calc_layers' function)
INPUT_FN = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'optimization/v0.4.0/flights/flight_20161109/interim_data/'...
    'SEAT2010_4.mat'];

% Path for output of saved manully traced layers
OUTPUT_FN = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'optimization/v0.4.0/flights/flight_20161109/man_layers/'...
    'SEAT2010_4.mat'];

% Load radar to workspace
radar = load(INPUT_FN);

% Run drawing routine
[man_layers] = draw_manual(radar);

%% Save manual layer output to disk for later use

% Save processed radar structure for future use
save(OUTPUT_FN, 'man_layers')