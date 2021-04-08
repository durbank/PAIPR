% Re-usable script to manually draw layers from radar echograms (using the
% 'draw_manual' function)

%% Set up environment

% Define data directory for flight path of interest
DATA_DIR = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'optimization/v0.4.0/flights/flight_20161109/'];

% Define name of site of interest
fn = 'Site_C.mat';



% Import src functions and dependencies
ROOT_DIR = fileparts(fileparts(pwd));
addpath(genpath(fullfile(ROOT_DIR, 'src/PAIPR-core/')))
addpath(genpath(fullfile(ROOT_DIR, 'src/Dependencies/')))
% addpath(genpath(fullfile(ROOT_DIR, 'src/optimize/')))

%% Initialize manual drawing function

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