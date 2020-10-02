% Script to optimize logistic regression coefficients for layer-likelihood
% calculations within the PAIPR algorithms

%% Add dependencies

ROOT_DIR = fileparts(fileparts(pwd));
DATA_DIR = ['/media/durbank/WARP/Research/Antarctica/Data/IceBridge/'...
    'optimization/v0.4.0/'];
NSIM = 225;

% Add PAIPR-core functions to path
addon_struct = dir(fullfile(ROOT_DIR, 'src/', 'PAIPR-core*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add PAIPR optimization functions to path
addon_struct = dir(fullfile(ROOT_DIR, 'src/', 'optimize*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile(ROOT_DIR, 'src/Dependencies/', ...
    'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

%%

% Get list of flightline directories
flight_list = dir(fullfile(DATA_DIR, 'flights'));
flight_list = flight_list(~ismember({flight_list.name},{'.','..'}));

for i=1:length(flight_list)
    
    content = dir(fullfile(flight_list(i).folder, flight_list(i).name, ...
        'interim_data'));
    f_dir = [content.isdir];
    sites = content(~f_dir);
    
    for j=1:length(sites)
        
        % Assign output path for results
        out_path = fullfile(DATA_DIR, 'params', ...
            strcat('params_', flight_list(i).name, '_', sites(j).name));
        
        % Load current data
        radar = load(fullfile(sites(j).folder, sites(j).name));
        tmp = load(fullfile(fileparts(sites(j).folder), 'man_layers', ...
            sites(j).name));
        man_layers = tmp.man_layers;
        
        % Create logical matrix of manual layers
        man_idx = cellfun(@(x) sub2ind(size(radar.data_smooth), ...
            round(x(:,2)),x(:,1)), man_layers, 'UniformOutput', false);
        man_grid = zeros(size(radar.data_smooth));
        for k = 1:length(man_idx)
            man_grid(man_idx{k}) = 1;
        end
        
        % Extract variables (to avoid unnecessary broadcasting)
        DB = radar.DB;
        depth = radar.depth;
        
        % Preallocate arrays
        r_params = zeros(1, size(DB,2));
        k_params = zeros(1, size(DB,2));
        SSE = zeros(1, size(DB,2));
        
        parfor k = 1:length(r_params)
            
            [r_params(k), k_params(k), SSE(k)] = opt_param(man_grid(:,k), ...
                DB(:,k), depth, 3, -8);
            
        end
        
        % Calculate logistic parameters using SSE-weighted sums
        weights = (1./SSE)./sum(1./SSE);
        r = sum(weights.*r_params);
        k = sum(weights.*k_params);
        
        % Calculate radar ages
        [radar, N_new] = radar_age(radar, r, std(r_params), k, NSIM);
        
        % Preallocate matrix of manual layer ages
        man_ages = zeros(size(radar.data_smooth));
        
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
        
        % Save the parameter distributions
        save(out_path, 'r_params', 'k_params', 'SSE')
        
        % Output path for final radar data
        radar_out = fullfile(DATA_DIR, 'flights', flight_list(i).name, ...
            'final_data', sites(j).name);
        radar.man_grid = man_grid;
        radar.man_ages = man_ages;
        save(radar_out, '-struct', 'radar')
        
    end
    
end
