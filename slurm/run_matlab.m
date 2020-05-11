

% Add PAIPR-core functions to path
addon_struct = dir(fullfile('src/', 'PAIPR-core*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add cresis-matlab-reader functions to path
addon_struct = dir(fullfile('src/dependencies/', ...
    'cresis-L1B-matlab-reader*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile('src/dependencies/', ...
    'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

%%

% If required, start parellel pool
poolobj=parpool('local',8);

% Run script
[success_codes] = process_SLURM(DATADIR, RHO_PATH, OUTDIR);
% Save sucess codes as .csv table

% If required, close parellel pool
delete(poolobj)

% Exit matlab
exit
