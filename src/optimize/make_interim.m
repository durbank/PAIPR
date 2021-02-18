% Script to generate interim data for drawing manual layers on processed
% echogram data

%% Add dependencies

ROOT_DIR = fileparts(fileparts(pwd));
addpath(genpath(fullfile(ROOT_DIR, 'src/PAIPR-core/')))
addpath(genpath(fullfile(ROOT_DIR, 'src/Dependencies/')))
addpath(genpath(fullfile(ROOT_DIR, 'src/optimize/')))

%% Generate interim data

% Assign data directory and density file
DATA_DIR = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'optimization/v0.4.0/flights/flight_20111109/'];
RHO_FILE = ['/media/durbank/WARP/Research/Antarctica/Data/CHPC/'...
    'flight-density/rho_20111109.csv'];

% Import density profiles and format as nested table
rho_table = import_rho(RHO_FILE);

% Get list of raw data site names
sites = dir(fullfile(DATA_DIR, 'raw_data'));
sites(ismember( {sites.name}, {'.', '..'})) = [];
site_names = {sites.name};

% Get list of existing interim datasets
datum = dir(fullfile(DATA_DIR, 'interim_data'));
datum(ismember( {datum.name}, {'.', '..'})) = [];
data_names = {datum.name};

for i=1:length(sites)
    
    % Create site file name
    site_fn = strcat(site_names{i}, '.mat');
    
    % Check if interim data already exist
    if ~ismember(site_fn, data_names)
        
        % Select raw data directory for current iteration
        raw_data = fullfile(DATA_DIR, 'raw_data', site_names{i});
        
        % PAIPR processing
        [radar] = manual_prep(raw_data, rho_table);

        % Export to relevant interim data directory
        output_path = fullfile(DATA_DIR, 'interim_data', site_fn);
        save(output_path, 'radar')
        
    else
        sprintf(...
            'Interim data already exists for %s...moving to next site', ...
            site_names{i})
    end
end







