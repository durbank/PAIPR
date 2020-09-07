% Script to optimize logistic regression coefficients for layer-likelihood
% calculations within the PAIPR algorithms


% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'G:/Research/Antarctica/Data/';
        addon_path = 'N:/MATLAB/Add-ons/';
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile(addon_path, 'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add export_fig to path
addon_struct = dir(fullfile(addon_path, 'altmany-export_fig*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add CReSIS OIB MATLAB reader functions to path
addon_struct = dir(fullfile(addon_path, 'cresis-L1B-matlab-readers*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add PAIPR-core functions to path
PAIPR_path = fullfile(cd,'..', 'PAIPR-core');
addpath(genpath(PAIPR_path))

% Add PAIPR optimization functions to path
opt_path = fullfile(cd, '..', 'optimize');
addpath(genpath(opt_path))

% Select path to directoryin containing training datasets
path_new = fullfile(data_path, "PAIPR-results/v0.3.0/optim");
% path_new = uigetdir(data_path, ...
%     'Select optimization directory to use');


%% Process raw OIB echograms with PAIPR for use in log regression optimization


% Get list of test core sites
f_list = dir(fullfile(path_new, 'raw_data'));
f_dir = [f_list.isdir];
f_named = ~ismember({f_list(:).name},{'.','..'});
sites = f_list(f_dir & f_named);


% Define number of Monte Carlo simulations to perform
Ndraw = 100;


% Select modeled density subset .mat file to use, and load to workspace
rho_fn = "rho_20111109subset.mat";
rho_path = "/media/durbank/WARP/Research/Antarctica/Data/PAIPR-results/v0.3.0/";
% [rho_fn, rho_path] = uigetfile([data_path '*.mat'], ...
%     "Select modeled density .mat subset");
rho_subset = load(fullfile(rho_path, rho_fn));
rho_subset = rho_subset.rho_subset;


for h=1:length(sites)
    
    % Define core site name
    core_site = sites(h).name;
    
    % Define path to interim PAIPR-processed data
    f_interim = fullfile(path_new, "interim_data", ...
        strcat(core_site, ".mat"));
    
    % Only process echograms without previous interim processing
    if ~isfile(f_interim)
        
        % Define directory containing raw echogram files
        echo_dir = fullfile(sites(h).folder, sites(h).name);
        
        
        % Get grouping indices of echograms
        [files, start_idx, end_idx] = OIB_chunk(echo_dir);
        
        
        
        
        
        
        horz_res = 25;
        
        
        echo_0 = orderfields(import_radar(...
            OIB_convert(fullfile(files(1).folder, files(1).name))));
        fields = fieldnames(echo_0);
        for i = 1:length(fields)
            echo_0.(fields{i}) = [];
        end
        
        
        
        %         for i=1:length(end_idx)
        i=1;
        
        echo_i = echo_0;
        for j = start_idx(i):end_idx(i)
            
            OIB_j = OIB_convert(fullfile(files(j).folder, files(j).name));
            struct_j = orderfields(import_radar(OIB_j));
            
            echo_i = cell2struct(cellfun(@horzcat, struct2cell(echo_i), ...
                struct2cell(struct_j), 'uni', 0), fieldnames(echo_i), 1);
            
            
        end
        echo_i.dist = pathdistps(echo_i.lat, echo_i.lon);
        
        %%
        
        [radar] = radar_stack(echo_i, horz_res);
        
        
        % Load modeled depth-density data from stats model output at specified
        % Easting/Northing coordinates
        [rho_data] = load_rho(rho_subset, radar.Easting, radar.Northing);
        
        % Convert to depth
        [radar] = radar_depth(radar, rho_data);
        
        % Calculate radar age-depth profile distributions (includes processing
        %signal-noise, radon transforms, layer tracing, likelihood assignments,
        % and age calculations)
        [radar] = calc_layers(radar, 'stream');
        
        %%
        
        
        man_file = fullfile(path_new, 'man_layers', ...
            strcat(core_site, '.mat'));
        if ~isfile(man_file)
            
            % Call function to trace manual layers and save as .mat file
            [man_layers] = draw_manual(radar, man_file);
        else
            
            man_layers = load(fullfile(path_new, 'man_layers', core_site));
            man_layers = man_layers.man_layers;
        end
        
        
        %% Tmp statements to deal with mismatched manually-traced layers
        
        if strcmp(core_site, 'SEAT2010_5')
            for k = 1:length(man_layers)
                vals = man_layers{k};
                vals(:,1) = vals(:,1)+200;
                man_layers{k} = vals;
            end
        elseif strcmp(core_site, 'SEAT2010_6')
            for k = 1:length(man_layers)
                vals = man_layers{k};
                vals(:,1) = vals(:,1)+185;
                man_layers{k} = vals;
            end
        end
        %%
        
        % Find the
        depth_max = ceil(max(cellfun(@(x) max(x(:,2)), man_layers)));
        
        man_idx = cellfun(@(x) ...
            sub2ind([depth_max, size(radar.data_smooth,2)], ...
            round(x(:,2)),x(:,1)), man_layers, 'UniformOutput', false);
        man_grid = zeros([depth_max size(radar.data_smooth,2)]);
        for k = 1:length(man_idx)
            man_grid(man_idx{k}) = 1;
        end
        
        
        man_Xbnds = [floor(min(cellfun(@(x) min(x(:,1)), man_layers))) ...
            ceil(max(cellfun(@(x) max(x(:,1)), man_layers)))];
        Xstart = max([1 man_Xbnds(1)]);
        Xend = min([size(radar.data_smooth,2) man_Xbnds(2)]);
        depth_end = min([depth_max length(radar.depth)]);
        
        radar.collect_time = radar.collect_time(Xstart:Xend);
        radar.Easting = radar.Easting(Xstart:Xend);
        radar.Northing = radar.Northing(Xstart:Xend);
        radar.Lat = radar.Lat(Xstart:Xend);
        radar.Lon = radar.Lon(Xstart:Xend);
        radar.elev = radar.elev(Xstart:Xend);
        radar.dist = radar.dist(Xstart:Xend);
        radar.dist = radar.dist-radar.dist(1);
        radar.depth = radar.depth(1:depth_end);
        radar.data_smooth = radar.data_smooth(1:depth_end,Xstart:Xend);
        radar.IM_grad = radar.IM_grad(1:depth_end,Xstart:Xend);
        radar.IM_QC = radar.IM_QC(1:depth_end,Xstart:Xend);
        radar.DB = radar.DB(1:depth_end,Xstart:Xend);
        radar.man_grid = man_grid(1:depth_end,Xstart:Xend);
        
        man_layers = load(fullfile(path_new, 'man_layers', core_site));
        man_layers = man_layers.man_layers;
        radar.man_layers = man_layers;
        
        save(fullfile(path_new, "interim_data", strcat(core_site, ".mat")), ...
            '-struct', 'radar')
    end
    
    %     end
    
end


%%
clearvars -except data_path path_new

% Get list of test core sites
f_list = dir(fullfile(path_new, 'interim_data'));
f_dir = [f_list.isdir];
sites = f_list(~f_dir);

for i = 1:length(sites)
    
    core_site = sites(i).name;
    output_path = fullfile(path_new, strcat("params_", core_site));
    
    if ~isfile(output_path)
        
        radar = load(fullfile(path_new, "interim_data", core_site));
        
        DB = radar.DB;
        depth = radar.depth;
        man_grid = radar.man_grid;
        
        r_params = zeros(1, size(DB,2));
        k_params = zeros(1, size(DB,2));
        SSE = zeros(1, size(DB,2));
        
        parfor j = 1:length(r_params)
            
            [r_params(j), k_params(j), SSE(j)] = opt_param(man_grid(:,j), ...
                DB(:,j), depth, 3, -8);
            
        end
        
        
        
        
        %         % Calculate logistic parameters using SSE-weighted sums
        %         weights = (1./SSE)./sum(1./SSE);
        %         r = sum(weights.*r_params);
        %         k = sum(weights.*k_params);
        %
        % Find radar ages
        Ndraw = 100;
        r = -6.8;
        k = 3;
        [radar] = radar_age(radar, r, k, Ndraw);
        
        
        % Preallocate matrix of manual layer ages
        man_ages = zeros(size(radar.data_smooth));
        
        % Define surface age and the year associated with the first pick of the
        % algorithm
        yr_vec = datevec(radar.collect_time(i));
        yr_pick1 = yr_vec(1);
        age_top = yr_vec(1) + (30*yr_vec(2)+yr_vec(3))/365;
        
        % Interpolate manual layer age-depth scales for each trace
        for k = 1:size(man_ages,2)
            depths_k = [0; radar.depth(logical(radar.man_grid(:,k)))];
            yrs_k = ([age_top yr_pick1:-1:yr_pick1-length(depths_k)+2])';
            man_ages(:,k) = interp1(depths_k, yrs_k, radar.depth, ...
                'linear', 'extrap');
        end
        
        % Diagnostic plot of manual and PAIPR ages
        plt_idx = randi(size(radar.data_smooth, 2));
        plt_idx = 500;
        
        figure
        hold on
        plot(radar.depth, median(squeeze(radar.ages(:,plt_idx,:)),2), 'r')
        plot(radar.depth, median(squeeze(radar.ages(:,plt_idx,:)),2) ...
            + std(squeeze(radar.ages(:,plt_idx,:)), [],2), 'r--')
        plot(radar.depth, median(squeeze(radar.ages(:,plt_idx,:)),2) ...
            - std(squeeze(radar.ages(:,plt_idx,:)), [],2), 'r--')
        plot(radar.depth, man_ages(:,plt_idx), 'b')
        hold off
        
        
        
        %%% Save the parameter distributions
        
        save(output_path, 'r_params', 'k_params', 'SSE')
        
    end
    
end