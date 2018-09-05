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


%% Comparison of trends across divide

% run import functions from accum_results

smb_oib = cellfun(@(x) median(median(x)), radar_OIB.SMB);
p = cellfun(@(x,y) robustfit(x, median(y,2)), radar_OIB.SMB_yr, radar_OIB.SMB, ...
    'UniformOutput', false);
trend_oib = cellfun(@(x) x(2), p);
change_oib = trend_oib./smb_oib;



inputs = {'SEAT10_4', 'SEAT10_5', 'SEAT10_6'};

change_seat = zeros(1, length(inputs));
seat_Easting = zeros(size(change_seat));
seat_elev = cores.elev(3:5);
for i = 1:length(inputs)
    name = inputs{i};
    file = fullfile(data_path, 'radar/SEAT_Traverses/results_data', ...
        strcat('grid', name, '.mat'));
    radar_i = load(file);
    smb1 = cellfun(@(x) median(median(x)), radar_i.SMB);
    p = cellfun(@(x,y) robustfit(x, median(y,2)), radar_i.SMB_yr, radar_i.SMB, ...
        'UniformOutput', false);
    trend1 = cellfun(@(x) x(2), p);
    change_seat(i) = median(trend1./smb1);
    seat_Easting(i) = mean(radar_i.Easting);
end

figure
hold on
plot(radar_OIB.Easting, change_oib, 'm.', 'MarkerSize', 10)
plot(seat_Easting, change_seat, 'rx', 'MarkerSize', 25)
hold off

figure
hold on
plot(seat_elev, change_seat, 'ro')
% label(seat_elev, change_seat, inputs)


