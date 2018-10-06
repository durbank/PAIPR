figure
imagesc(peaks_raw)
hold on
for k = 1:length(layers)
    [r,c] = ind2sub(size(peaks_raw), layers{k});
    plot(c, r, '.', 'MarkerSize', 15)
end
hold off

figure
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
hold on
for k = 1:length(layers_idx)
    [r,c] = ind2sub(size(radar.data_smooth), layers_idx{k});
    r_scale = r*0.02;
    c_scale = c*mean(diff(radar.dist));
    plot(c_scale, r_scale, 'LineWidth', 3)
end
hold off

%% Map concat (combine N-E points of multiple files to determine map position

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

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

path2 = 'radar/SEAT_Traverses/SEAT2010Kuband/ProcessedSEAT2010/grid_SEAT10_6';
radar_dir = fullfile(data_path, path2);

% Get list of radar files in directory
wild = '*.mat';
files = dir(fullfile(radar_dir, wild));

figure
hold on
h1 = scatter(cores.Easting, cores.Northing, 100, 'b', 'filled');
h2 = plot(radar.Easting(1), radar.Northing(1), 'r', 'LineWidth', 2);     % Correctly display radar as line in legend
plot(radar.Easting, radar.Northing, 'r.', 'MarkerSize', 0.10)

for i = 1:length(files)
    data = radar_clean(fullfile(radar_dir, files(i).name));
    plot(data.Easting, data.Northing, '.')
end
hold off


%% Compare std of trace distributiions with those nearby

cores_loop = {'SEAT10_4' 'SEAT10_5' 'SEAT10_6'};
near_dist = 6500;
i = 3;
core_i = cores.(cores_loop{i});

D_SEATi = pdist2([core_i.Easting, core_i.Northing], [SEAT_E' SEAT_N']);
SEATi_idx = D_SEATi <= near_dist;
[~, SEATi_near] = min(D_SEATi(SEATi_idx));

SEAT_data = SEAT_SMB_MC(SEATi_idx);
SMB_length = min(cellfun(@(x) size(x,1), SEAT_data));
SEAT_std = cell2mat(cellfun(@(x) std(x(1:SMB_length,:),[],2), ...
    SEAT_data, 'UniformOutput', 0));
figure
yr_start = SEAT_yr{1}(1);
plot(yr_start:-1:(yr_start-SMB_length+1), SEAT_std, 'LineWidth', 0.5)
hold on
plot(yr_start:-1:(yr_start-SMB_length+1), SEAT_std(:,SEATi_near), ...
    'r', 'LineWidth', 3)
hold off

D_OIBi = pdist2([core_i.Easting, core_i.Northing], [OIB_E' OIB_N']);
OIBi_idx = D_OIBi <= near_dist;
[~, OIBi_near] = min(D_OIBi(OIBi_idx));

OIB_data = OIB_SMB_MC(OIBi_idx);
SMB_length = min(cellfun(@(x) size(x,1), OIB_data));
OIB_std = cell2mat(cellfun(@(x) std(x(1:SMB_length,:),[],2), ...
    OIB_data, 'UniformOutput', 0));
figure
yr_start = OIB_yr{1}(1);
plot(yr_start:-1:(yr_start-SMB_length+1), OIB_std, 'LineWidth', 0.5)
hold on
plot(yr_start:-1:(yr_start-SMB_length+1), OIB_std(:,OIBi_near), ...
    'r', 'LineWidth', 3)
hold off

