figure
imagesc(peaks_raw)
hold on
for k = 1:length(layers)
    [r,c] = ind2sub(size(peaks_raw), layers{k});
    plot(c, r, '.', 'MarkerSize', 15)
end
hold off

figure
imagesc(radar.data_smooth, [-2 2])
hold on
for k = 1:length(radar.layers)
    [r,c] = ind2sub(size(radar.data_smooth), radar.layers{k});
    plot(c, r, 'LineWidth', 3)
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

%%

Ndraw = seg1.Ndraw;
cores = seg1.cores;
radar = struct('collect_date', seg1.radar.collect_date, 'Easting', ...
    [seg1.radar.Easting(100:end-100) seg2.radar.Easting(100:end-100) seg3.radar.Easting(100:end-100)],...
    'Northing', [seg1.radar.Northing(100:end-100) seg2.radar.Northing(100:end-100) seg3.radar.Northing(100:end-100)],...
    'dist', [seg1.radar.dist(100:end-100) seg2.radar.dist(100:end-100) seg3.radar.dist(100:end-100)],...
    'depth', seg1.radar.depth, ...
    'data_smooth', [seg1.radar.data_smooth(:,100:end-100) seg2.radar.data_smooth(:,100:end-100) seg3.radar.data_smooth(:,100:end-100)],...
    'likelihood', [seg1.radar.likelihood(:,100:end-100) seg2.radar.likelihood(:,100:end-100) seg3.radar.likelihood(:,100:end-100)],...
    'ages', [seg1.radar.ages(:,100:end-100,:) seg2.radar.ages(:,100:end-100,:) seg3.radar.ages(:,100:end-100,:)]);
radar.SMB_yr = [seg1.radar.SMB_yr(100:end-100) seg2.radar.SMB_yr(100:end-100) seg3.radar.SMB_yr(100:end-100)];
radar.SMB = [seg1.radar.SMB(100:end-100) seg2.radar.SMB(100:end-100) seg3.radar.SMB(100:end-100)];



%% 
% Trace idx to investigate
i = randi(size(radar.data_smooth, 2));

% Find the nearest cores to the radar data (for comparison plots)
[~, cores_near_idx] = sort(pdist2([radar.Easting(i) radar.Northing(i)], ...
    [cores.Easting' cores.Northing'], 'Euclidean'));
core_near1 = cores.(cores.name{cores_near_idx(1)});
core_near2 = cores.(cores.name{cores_near_idx(2)});
core_near3 = cores.(cores.name{cores_near_idx(3)});

% Calculate the mean age-depth scale and std for radar trace i
age_mean = mean(squeeze(radar.ages(:,i,:)), 2);
age_std = std(squeeze(radar.ages(:,i,:)), [], 2);

% Plot full radargram
yr_idx = logical([diff(floor(age_mean)); 0]);
depth = radar.depth(yr_idx);
col = i*ones(length(depth),1);
% row = find(radar.likelihood(:,i)>0.75);
% col = i*ones(length(row),1);
figure('Position', [200 200 1500 800])
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
colorbar
xlabel('Distance along profile (m)')
ylabel('Depth (m)')
hold on
plot(radar.dist(col), depth, 'r.', 'MarkerSize', 25)
xlim([0 radar.dist(end)])
ylim([0 radar.depth(end)])
set(gca, 'Ydir', 'reverse', 'FontSize', 18)
hold off

% Age-depth scale comparison between radar trace and nearest cores
figure
hold on
h1 = plot(core_near1.depth, mean(core_near1.ages, 2), 'b', 'LineWidth', 2);
plot(core_near1.depth, mean(core_near1.ages, 2) + 2*std(core_near1.ages, [], 2), 'b--')
plot(core_near1.depth, mean(core_near1.ages, 2) - 2*std(core_near1.ages, [], 2), 'b--')
h2 = plot(core_near2.depth, core_near2.age, 'c', 'LineWidth', 2);
h3 = plot(core_near3.depth, core_near3.age, 'c--', 'LineWidth', 1);
h4 = plot(radar.depth, age_mean, 'r', 'LineWidth', 2);
plot(radar.depth, age_mean + 2*age_std, 'r--', 'LineWidth', 0.5)
plot(radar.depth, age_mean - 2*age_std, 'r--', 'LineWidth', 0.5)
ylabel('Calendar Year')
xlabel('Depth (m)')
legend([h1 h2 h3 h4], 'Nearest core age (manual)', '2nd nearest core', ...
    '3rd nearest core', 'Radar age (automated)', 'Location', 'ne')
set(gca, 'FontSize', 10)
hold off
