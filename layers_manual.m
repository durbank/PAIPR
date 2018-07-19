


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

%%

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
        depth_interp = (0:core_res:radar.depth(end,i))';
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

%%

num_layer = sum(logical(radar.man_layers));
mov_num = num_layer - movmedian(num_layer, 25);
rm_idx = mov_num < 0;

radar.Easting(rm_idx) = [];
radar.Northing(rm_idx) = [];
radar.man_layers(:,rm_idx) = [];

% Define surface age and the year associated with the first pick of the 
% algorithm
age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);
yr_idx = logical(radar.man_layers);

ages = nan(size(yr_idx));
err_idx = [];
for i = 1:size(yr_idx, 2)
    
    depths_i = [0; radar.depth(yr_idx(:,i))];
    yrs_i = ([age_top yr_pick1:-1:yr_pick1-length(depths_i)+2])';
    try
        ages(:,i) = interp1(depths_i, yrs_i, radar.depth, 'linear');
%         ages(:,i) = interp1(depths_i, yrs_i, radar.depth, 'linear', 'extrap');
    catch
        sprintf('Error in age interpolation for trace %u. Removing trace from data', i)
        err_idx = [err_idx i];
    end
end


radar.Easting(err_idx) = [];
radar.Northing(err_idx) = [];
radar.man_layers(:,err_idx) = [];
ages(:,err_idx) = [];
radar.ages = ages;

radar = rmfield(radar, {'dist', 'data_stack', 'rho_coeff', 'rho_var'});

fn = 'SEAT10_manual_layers.mat';
output_path = fullfile(data_path, 'radar/SEAT_Traverses/results_data', fn);
save(output_path, '-struct', 'radar', '-v7.3')