% Script to run PAIPR analyses locally (as opposed to HPC via SLURM batch
% scripts)

%% Add dependencies

ROOT_DIR = fileparts(fileparts(pwd));

% Add PAIPR-core functions to path
addon_struct = dir(fullfile(ROOT_DIR, 'src/', 'PAIPR-core*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add cresis-matlab-reader functions to path
addon_struct = dir(fullfile(ROOT_DIR, 'src/Dependencies/', ...
    'cresis-L1B-matlab-reader*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile(ROOT_DIR, 'src/Dependencies/', ...
    'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

%% 

DATADIR = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'manual_layers/SEAT10_4/raw_data/'];
RHOFILE = ['/media/durbank/WARP/Research/Antarctica/Data/CHPC/'...
    'flight-density/rho_20111109.csv'];
OUTDIR = '/home/durbank/scratch/';
NSIM = 500;

% [success_codes] = process_PAIPR(DATADIR, RHOFILE, OUTDIR, NSIM);
[success_codes] = process_PAIPR(...
    DATADIR, RHOFILE, OUTDIR, NSIM, 'VerboseOutput', true);