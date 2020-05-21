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
    'WAIS-central/20111109/'];
RHOFILE = ['/media/durbank/WARP/Research/Antarctica/Data/rho-stats/Data/'...
    'cent_interior/Data/densities/rho_20111109.csv'];
OUTDIR = '/home/durbank/Documents/scratch/';
NSIM = 100;

[success_codes] = process_SLURM(DATADIR, RHOFILE, OUTDIR, NSIM);
% [success_codes] = process_SLURM(...
%     DATADIR, RHOFILE, OUTDIR, NSIM, 'VerboseOutput', true);