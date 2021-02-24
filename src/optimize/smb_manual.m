% Script to generate and save SMB estimates based on manually traced layers

%% Set environment

ROOT_DIR = fileparts(fileparts(pwd));
DATA_DIR = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'optimization/v0.4.0/'];
RHOFILE = ['/media/durbank/WARP/Research/Antarctica/Data/CHPC/'...
    'flight-density/rho_20111109.csv'];
OUT_DIR = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'optimization/v0.4.0/SMB_results/'];
NSIM = 100;

% Add dependencies and functions to path
addpath(genpath(fullfile(ROOT_DIR, 'src')))

%% Generate file names to loop through

% Get list of flightline directories
flight_list = dir(fullfile(DATA_DIR, 'flights'));
flight_list = flight_list(~ismember({flight_list.name},{'.','..'}));

% Import density profiles and format as nested table
rho_table = import_rho(RHOFILE);

%%

for i=1:length(flight_list)
    
    % Get name of flight
    flight_parts = strsplit(flight_list(i).name, '_');
    flight_name = flight_parts{end};
    
    % Get list of interim file paths
    content = dir(fullfile(flight_list(i).folder, flight_list(i).name, ...
        'interim_data'));
    f_dir = [content.isdir];
    files_interim = content(~f_dir);
    
    % Get list of manual layer files
    content = dir(fullfile(flight_list(i).folder, flight_list(i).name, ...
        'man_layers'));
    f_dir = [content.isdir];
    files_manual = content(~f_dir);
    
    parfor j=1:length(files_interim)
        
        % Load current radar and man_layer files
        radar = load(fullfile(...
            files_interim(j).folder, files_interim(j).name));
        tmp = load(fullfile(files_manual(j).folder, files_manual(j).name));
        man_layers = tmp.man_layers;
        
        % Ensure proper overlap between radar and man_layers
        [radar] = optim_format(radar, man_layers);
        
        % Construct updated matrix of manual layers based on man_layer 
        % indices
        man_grid = zeros(size(radar.data_smooth));
        for k = 1:length(radar.man_layers)
            man_grid(radar.man_layers{k}) = k;
        end
        
        % Preallocate matrix of manual layer ages
        man_ages = zeros(size(man_grid));
        
        % Define surface age and the year associated with the first pick
        % of the algorithm
        yr_vec = datevec(mean(radar.collect_time));
        yr_pick1 = yr_vec(1);
        age_top = yr_vec(1) + (30*yr_vec(2)+yr_vec(3))/365;
        
        % Interpolate manual layer age-depth scales for each trace
        for k = 1:size(man_ages,2)
            depths_k = [0; radar.depth(logical(man_grid(:,k)))];
            yrs_k = ([age_top yr_pick1:-1:yr_pick1-length(depths_k)+2])';
            man_ages(:,k) = interp1(depths_k, yrs_k, radar.depth, ...
                'linear', 'extrap');
        end
        
        % Reshape ages to match Ndraw and add to radar structure
        radar.ages = repmat(man_ages, 1, 1, NSIM);
        
        % Load modeled depth-density data from stats model output at
        % specified Easting/Northing coordinates
        [rho_data] = load_rho(...
            rho_table, radar.Easting, radar.Northing);
        
        % Calculate annual accumulation
        [radar] = calc_SWE(radar, rho_data, NSIM);
        
        % Perform QC check on echogram image
        [QC_med, QC_val, QC_flag, depth_idx, yr_cutoff] = QC_check(...
            radar, 0.50, 0.10);
        radar.QC_med = QC_med;
        radar.QC_val = QC_val;
        radar.QC_flag = QC_flag;
        radar.QC_depth_idx = depth_idx;
        radar.QC_yr = yr_cutoff;
        
        % Generate file names and paths under which to save data
        [~, filename, ~] = fileparts(fullfile(...
            files_interim(j).folder, files_interim(j).name));
        csv_output = fullfile(OUT_DIR, flight_name, ...
            strcat(filename, '_', flight_name, '.csv'));
        
        % Save csv outputs to disk
        [success_code] = parsave(radar, csv_output, 'mixture', 200);
        
    end
    
end
