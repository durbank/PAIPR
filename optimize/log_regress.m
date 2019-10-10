% Script to estimate the logistic regression coefficients for calculating
% layer likelihood scores

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/Data/';
        addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
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


%% Process raw OIB echograms with PAIPR for use in log regression optimization

% % Define number of Monte Carlo simulations to perform
% Ndraw = 100;
%
% % Load core data from file (data used was previously generated using
% % import_cores.m)
% core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
% cores = load(core_file);
%
% % Select directory containing raw OIB echograms to process
% input_dir = uigetdir(fullfile(data_path, "PAIPR-results/vWIP", ...
%     "coeffs_branch/manual_layers"), ...
%     "Select directory containing input files");
%
% % Select output directory in which to save processed echogram and manual
% % layers
% [output_dir] = uigetdir(input_dir, ...
%     "Select directory to output processed echogram");
%
% % Process OIB echogram with PAIPR
% [radar] = PAIPR_draw(input_dir, cores, Ndraw);
%
% %%%
%
% % Clip depth-related variables to final cutoff depth
% cutoff = 25;
% core_res = 0.02;
% cut_idx = min([(round(cutoff/core_res)+1) length(radar.depth)]);
% radar_tmp = struct('collect_date', radar.collect_date, 'Easting', ...
%     radar.Easting, 'Northing', radar.Northing, 'dist', radar.dist, ...
%     'depth', radar.depth(1:cut_idx), ...
%     'data_smooth', radar.data_smooth(1:cut_idx,:),...
%     'peaks', radar.peaks(1:cut_idx,:), 'groups', radar.groups(1:cut_idx,:));
%
% % Clip 5 km off the start/end of the processed radargram (in order to
% % properly match location and dimensions of manually traced layers
% overlap = 10000;
% horz_res = 25;
% clip = round(0.5*overlap/horz_res);
% radar = struct('collect_date', radar_tmp.collect_date, ...
%     'Easting', radar_tmp.Easting(clip:end-clip),...
%     'Northing', radar_tmp.Northing(clip:end-clip), ...
%     'dist', radar_tmp.dist(clip:end-clip), 'depth', radar_tmp.depth, ...
%     'data_smooth', radar_tmp.data_smooth(:,clip:end-clip),...
%     'peaks', radar_tmp.peaks(:,clip:end-clip), ...
%     'groups', radar_tmp.groups(:,clip:end-clip));
%
% % Redefine radargram distances based on clipped data
% radar.dist = radar.dist - radar.dist(1);
%
% % Find layer member indices based on new clipped record
% layers = cell(1, max(radar.groups(:)));
% for j = 1:length(layers)
%     layers{j} = find(radar.groups == j);
% end
% radar.layers = layers(~cellfun(@isempty, layers));
%
% %%%
%
% % Save processed radar structure for future use
% output_path = fullfile(output_dir, "PAIPR_out.mat");
% save(output_path, '-struct', 'radar', '-v7.3')


%%

% Select directory containing PAIPR-results and manual layers (picked using
% the `draw_manual.m` file)
input_dir = uigetdir(fullfile(data_path, "PAIPR-results/vWIP", ...
    "coeffs_branch/manual_layers"), ...
    "Select directory containing input files");

% Load PAIPR radar data
radar = load(fullfile(input_dir, "PAIPR_out"));

% Load manually picked layer indices for above radar echogram
load(fullfile(input_dir, "manual_layers"))

%% Calculate distance brightnesses

% Weighting coefficients for mean brightness and echogram distance
peak_w = 1/mean(radar.peaks(radar.peaks>0));
dist_w = 1/radar.dist(end);
peak_w = 1;
dist_w = 1;

horz_res = round(mean(diff(radar.dist)));
lay_dist = horz_res*cellfun(@length, radar.layers);

% Map distance-brightness values to the location within the radar
% matrix of the ith layer
DB = zeros(size(radar.peaks));
for i = 1:length(radar.layers)
    DB(radar.layers{i}) = peak_w*dist_w*...
        radar.peaks(radar.layers{i}).*lay_dist(i);
end

%%

man_idx = cellfun(@(x) sub2ind(size(radar.data_smooth), round(x(:,2)),x(:,1)), ...
    man_layers, 'UniformOutput', false);
man_peaks = zeros(size(radar.data_smooth));
for i = 1:length(man_idx)
    man_peaks(man_idx{i}) = 1;
end

%%

r_params = zeros(1, size(radar.data_smooth,2));
k_params = zeros(1, size(radar.data_smooth,2));
SSE = zeros(1, size(radar.data_smooth,2));
depth = radar.depth;

parfor i = 1:length(r_params)
    
    [r_params(i), k_params(i), SSE(i)] = opt_param(...
        man_peaks(:,i), DB(:,i), depth);
    
end


%%

output_path = fullfile(input_dir, "params_output.mat");
save(output_path, 'r_params', 'k_params', 'SSE')
