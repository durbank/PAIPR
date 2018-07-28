


% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'work';
%         computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
                
            case 'laptop'
                data_path = 'C:/Users/durba/Documents/Research/Antarctica/Data/';
                addon_path = 'C:/Users/durba/Documents/MATLAB/Addons/';
        end
        
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Number of Monte Carlo simulations
Ndraw = 100;

% Import firn core data
% [cores] = import_cores(fullfile(data_path, 'Ice-cores/SEAT_cores/', ...
%     'DGK_core_data.xlsx'), Ndraw);

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);


radar_file = fullfile(data_path, ...
    'radar/SEAT_Traverses/ManuallyTracedLayers/layers_ku_band_2010_manual.mat');

%% Modified from first portion of radar_age.m

% Conversion to depth
[radar] = radar_depth(radar_file, cores);

horz_res = 25;
[radar] = radar_stack(radar, horz_res);

% Stationarize the radar response by differencing traces with a smoothing 
% spline
s = zeros(size(radar.data_stack));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth(:,i), radar.data_stack(:,i), 0.95, radar.depth(:,i));
end
radar_stat = radar.data_stack - s;

% Remove linear trend in variance (attentuation with depth) and convert to
% standardized values (z-score statistics)
radar_Z = zeros(size(radar_stat));
for i = 1:size(radar_stat, 2)
    data = radar_stat(:,i);
    % Frame length to define local variance
    half_frame = round(0.5*length(data)/5); 
    var0 = movvar(data, 2*half_frame, 'EndPoints', 'discard');
    x = (half_frame:length(data)-half_frame)';
    % Linear trend in variance
    EQ = polyfit(x, var0, 1);
    x_mod = (1:length(data))';
    mod = polyval(EQ, x_mod);
    % Standardize variance-corrected data
    radar_Z(:,i) = data./sqrt(abs(mod));
end

% Define the vertical resolution of the core data
core_res = 0.02;

% Define the cutoff depth for radar traces and find index of crossover
% depth
cutoff = 25;
depth_bott = floor(min([min(radar.depth(end,:)) cutoff]));

% Trim radar traces to cutoff depth and interpolate data to vertical scale
% of the firn cores
radarZ_interp = zeros(depth_bott/core_res+1, size(radar.data_stack, 2));
for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:core_res:radar.depth(end,i));
    radarZ_i = interp1(radar.depth(:,i), radar_Z(:,i), depth_interp, 'pchip');
    radarZ_interp(:,i) = radarZ_i(1:size(radarZ_interp, 1));
end

% If manual layer picks are present, transform them to same depth and
% vertical scale as the interpolated radar data
if isfield(radar, 'man_layers')
    man_interp = zeros(size(radarZ_interp));
    for i = 1:size(radar.data_stack, 2)
%         depth_interp = (0:core_res:radar.depth(end,i))';
        layer_idx = logical(radar.man_layers(:,i));
        layer_num = radar.man_layers(layer_idx,i);
        man_depth_i = radar.depth(layer_idx,i);
        depth_idx = round(man_depth_i/core_res) + 1;
        man_interp(depth_idx(man_depth_i<=cutoff),i) = layer_num(man_depth_i<=cutoff);
    end
    radar.man_layers = man_interp;
end

% Assign structure output depth to interpolated depths
radar.depth = (0:core_res:depth_bott)';

% % Smooth the laterally averaged radar traces with depth based on a 3rd
% % order Savitzky-Golay filter with a window of 9 frames (~20 m)
% radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);

% Clear unnecessary variables
clearvars -except file cores Ndraw radar horz_res core_res data_path

%% Modified from final portion of radar_age.m

num_layer = sum(logical(radar.man_layers));
mov_num = num_layer - movmedian(num_layer, 25);
rm_idx = mov_num < 0;

radar.Easting(rm_idx) = [];
radar.Northing(rm_idx) = [];
radar.rho_coeff(:,rm_idx) = [];
radar.rho_var(:,rm_idx) = [];
radar.man_layers(:,rm_idx) = [];

% Define surface age and the year associated with the first pick of the 
% algorithm
% age_top = radar.collect_date;
% yr_pick1 = ceil(radar.collect_date - 1);
yr_pick1 = floor(radar.collect_date);
yr_idx = logical(radar.man_layers);

ages = nan(size(yr_idx));
err_idx = [];
for i = 1:size(yr_idx, 2)
    
%     depths_i = [0; radar.depth(yr_idx(:,i))];
%     yrs_i = ([age_top yr_pick1:-1:yr_pick1-length(depths_i)+2])';
    depths_i = radar.depth(yr_idx(:,i));
    yrs_i = (yr_pick1:-1:yr_pick1-length(depths_i)+1)';
    try
        ages(:,i) = interp1(depths_i, yrs_i, radar.depth, 'linear');
%         ages(:,i) = interp1(depths_i, yrs_i, radar.depth, 'linear', 'extrap');
    catch
        sprintf('Error in age interpolation for trace %u. Removing trace from data', i);
        err_idx = [err_idx i];
    end
end

radar.Easting(err_idx) = [];
radar.Northing(err_idx) = [];
radar.rho_coeff(:,err_idx) = [];
radar.rho_var(:,err_idx) = [];
radar.man_layers(:,err_idx) = [];
ages(:,err_idx) = [];
radar.ages = ages;

% %% Modified from calc_SWE.m
% 
% % Generate modelled std of density with depth for radar data
% rho_std = sqrt((radar.rho_var(1,:)-radar.rho_var(2,:))./...
%     (radar.rho_var(3,:).^radar.depth) + radar.rho_var(2,:));
% 
% % Generate density with depth model for radar data
% rho_mod = radar.rho_coeff(1,:).*radar.depth.^radar.rho_coeff(2,:) + ...
%     radar.rho_coeff(3,:);
% 
% %%% Calculate annual accumulation for each radar trace
% 
% % Calculate accumulation at each depth interval (0.02 m) with simulated
% % noise based on the variance in core density
% noise_rho = 1000*rho_std.*randn(size(radar.ages));
% accum_dt = 0.02*(1000*repmat(rho_mod, 1, 1, Ndraw) + noise_rho);
% 
% % Linearly interpolate traces with missing or problematic age-depth scales
% % sim_nan = all(isnan(radar.ages), 3);
% % trace_nan = any(sim_nan);
% radar.ages = fillmissing(radar.ages, 'linear', 2);
% 
% % Define age-depth profile based on the median of all MC age profiles and 
% % the std dev with depth of all MC age profiles for each trace (avoids  
% % integer year jumps in accumulation estimates and non-monotonically 
% % decreasing ages)
% age_std = squeeze(std(radar.ages, [], 3));
% age_noise = randn(1, 1, Ndraw).*age_std;
% ages = repmat(median(radar.ages, 3), 1, 1, Ndraw) + age_noise;
% % ages = radar.age;
% 
% 
% 
% % Define first year with complete accumulation data and earliest year with
% % observations within the data set
% yr_top = floor(max(max(ages(1,:,:))));
% yr_end = ceil(min(min(ages(end,:,:))));
% 
% % Define initial accumulation year vector (will be iteratively modified at
% % each trace in for loop)
% accum_yr_init = (yr_top-1:-1:yr_end)';
% 
% % Preallocation of cell arrays for accumulation years and annual
% % accumulation rate
% accum_yr = cell(1, size(radar.ages, 2));
% accum = cell(1, size(radar.ages, 2));
% for i = 1:size(accum, 2)
%     
%     % Preallocate vector for accumulation rate at ith trace
%     accum_i = zeros(length(accum_yr_init), Ndraw);
%     for j = 1:Ndraw
%         
%         % Calculate indices of integer ages for jth simulation of the ith
%         % trace (with added noise from uncertainty in exact point in time
%         % of the accumulation peak, using a std dev of 1 month)
%         years_j = ages(:,i,j);
%         yr_idx = logical([diff(floor(years_j)); 0]);
%         yr_loc = find(yr_idx);
%         loc_temp = yr_loc;
%         loc_temp(2:end-1) = yr_loc(2:end-1) + ...
%             round(movmean(1*diff(yr_loc(1:end-1))/12, 2).*...
%             randn(length(yr_loc)-2, 1));
%         loc_idx = loc_temp<1;
%         loc_temp(loc_idx) = yr_loc(loc_idx);
%         yr_loc = sort(loc_temp, 'ascend');
%         yr_loc(yr_loc>size(accum_dt, 1)) = size(accum_dt, 1);
%         
%         % Integrate accumulation at each depth point for each whole year in
%         % the jth simulation of the ith trace
%         n_length = length(yr_loc) - 1;
%         accum_j = zeros(n_length, 1);
%         for n = 1:n_length
%             accum_j(n) = sum(accum_dt(yr_loc(n)+1:yr_loc(n+1),i,j));
%         end
%         accum_i(1:n_length,j) = accum_j;
%     end
%     
%     % Output results to respective arrays
%     accum_idx = find(all(accum_i, 2), 1, 'last');
%     accum_clip = accum_i(1:accum_idx,:);
%     accum_yr{i} = accum_yr_init(1:accum_idx);
%     accum{i} = accum_clip;
% end
% 
% 
% %%% Output results to radar and core structures
% 
% radar.SMB_yr = accum_yr;
% radar.SMB = accum;

%% Final formatting and export to disk

radar = rmfield(radar, {'dist', 'data_stack', 'rho_coeff', 'rho_var'});

fn = 'SEAT10_manual_layers.mat';
output_path = fullfile(data_path, 'radar/SEAT_Traverses/results_data', fn);
save(output_path, '-struct', 'radar', '-v7.3')