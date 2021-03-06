% Script to compare the effect of differences in radargram length over the
% same area

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
addon_folder = fullfile(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))
% Add export_fig to path
addon_folder = fullfile(addon_path, 'altmany-export_fig-cafc7c5/');
addpath(genpath(addon_folder))

%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

Ndraw = 100;

% Load previously processed full OIB radar transect
% radar_full = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/OIB_SEAT10_4to10_6.mat'));


%% Concatenate adjacent OIB radar files (modified from OIB_concat.m)
%%% Using files 10-12 (this is approximately 15 km of data)
% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Assign directory of OIB data
OIB_dir = fullfile(data_path, 'IceBridge/SEAT10_4to10_6/2011_SNO');

% Get list of radar files in directory
wild = '*.nc';
files = dir(fullfile(OIB_dir, wild));

% Assign the data output directory (for saving processed files)
output_dir = fullfile('radar/SEAT_Traverses/results_data/comparisons/dataset_length');

%% radar_15km

% tic
% % Initialize for loop with first file in directory
% data_struct = OIB_import(fullfile(OIB_dir, files(10).name));
% fld_names = fieldnames(data_struct);
% data_cells = struct2cell(data_struct);
% 
% % Concatenate adjacent OIB radar files (modified from OIB_concat.m)
% for i = 11:12
%     data_i = struct2cell(OIB_import(fullfile(OIB_dir, files(i).name)));
%     data_cells = cellfun(@horzcat, data_cells, data_i, 'UniformOutput', 0);
% end
% 
% % Average values for TWTT, time_trace, and collect_date
% data_cells{8} = median(data_cells{8}, 2);
% data_cells{9} = median(data_cells{9});
% data_cells{10} = median(data_cells{10});
% 
% % Flip data (so that the processing direction matches SEAT direction from
% % SEAT10-4 to SEAT10-6)
% data_flip = cellfun(@fliplr, data_cells, 'UniformOutput', false);
% 
% radar_tmp = cell2struct(data_flip, fld_names, 1);
% radar_tmp.dist = pathdist(radar_tmp.lat, radar_tmp.lon);
% 
% toc
% % [radar_tmp] = OIB_age(radar_tmp, cores, Ndraw);
% [radar_tmp] = radar_RT(radar_tmp, cores, Ndraw);
% [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
% toc
% 
% %%% NEED TO FIX FIELD NAMES (NAME FIELDS EXPLICITLY TO AVOID
% %%% DIFFICULTLIES IN TRACING ERRORS IN THE FUTURE %%%
% clip = 100;
% fld_names = fieldnames(radar_tmp);
% radar_15km = struct(fld_names{1}, radar_tmp.(fld_names{1}), fld_names{2}, ...
%     radar_tmp.(fld_names{2})(clip:end-clip), ...
%     fld_names{3}, radar_tmp.(fld_names{3})(clip:end-clip),...
%     fld_names{4}, radar_tmp.(fld_names{4})(clip:end-clip),...
%     fld_names{5}, radar_tmp.(fld_names{5}), ...
%     fld_names{9}, radar_tmp.(fld_names{9})(:,clip:end-clip),...
%     fld_names{10}, radar_tmp.(fld_names{10})(:,clip:end-clip),...
%     fld_names{11}, radar_tmp.(fld_names{11})(:,clip:end-clip),...
%     fld_names{13}, radar_tmp.(fld_names{13})(:,clip:end-clip),...
%     fld_names{15}, radar_tmp.(fld_names{15})(:,clip:end-clip,:));
% 
% radar_15km.(fld_names{16}) = radar_tmp.(fld_names{16})(clip:end-clip);
% radar_15km.(fld_names{17}) = radar_tmp.(fld_names{17})(clip:end-clip);
% 
% % Save file in output directory
% flnm = 'OIB_15km.mat';
% output_path = fullfile(data_path, output_dir, flnm);
% save(output_path, '-struct', 'radar_15km', '-v7.3')
% 
% % Clear saved variables (reduces memory load on additional datasets)
% clearvars -except data_path OIB_dir files output_dir cores Ndraw
% toc
% Processing times: [6s 133s 143] 

radar_15km = load(fullfile(data_path, ...
    'radar/SEAT_Traverses/results_data/comparisons/dataset_length/OIB_15km.mat'));


%% radar_25km

% tic
% % Initialize for loop with first file in directory
% data_struct = OIB_import(fullfile(OIB_dir, files(9).name));
% fld_names = fieldnames(data_struct);
% data_cells = struct2cell(data_struct);
% 
% for i = 10:13
%     data_i = struct2cell(OIB_import(fullfile(OIB_dir, files(i).name)));
%     data_cells = cellfun(@horzcat, data_cells, data_i, 'UniformOutput', 0);
% end
% 
% % Average values for TWTT, time_trace, and collect_date
% data_cells{8} = median(data_cells{8}, 2);
% data_cells{9} = median(data_cells{9});
% data_cells{10} = median(data_cells{10});
% 
% % Flip data (so that the processing direction matches SEAT direction from
% % SEAT10-4 to SEAT10-6)
% data_flip = cellfun(@fliplr, data_cells, 'UniformOutput', false);
% 
% radar_tmp = cell2struct(data_flip, fld_names, 1);
% radar_tmp.dist = pathdist(radar_tmp.lat, radar_tmp.lon);
% toc
% 
% [radar_tmp] = radar_RT(radar_tmp, cores, Ndraw);
% [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
% toc
% 
% clip = 100;
% fld_names = fieldnames(radar_tmp);
% radar_25km = struct(fld_names{1}, radar_tmp.(fld_names{1}), fld_names{2}, ...
%     radar_tmp.(fld_names{2})(clip:end-clip), ...
%     fld_names{3}, radar_tmp.(fld_names{3})(clip:end-clip),...
%     fld_names{4}, radar_tmp.(fld_names{4})(clip:end-clip),...
%     fld_names{5}, radar_tmp.(fld_names{5}), ...
%     fld_names{9}, radar_tmp.(fld_names{9})(:,clip:end-clip),...
%     fld_names{10}, radar_tmp.(fld_names{10})(:,clip:end-clip),...
%     fld_names{11}, radar_tmp.(fld_names{11})(:,clip:end-clip),...
%     fld_names{13}, radar_tmp.(fld_names{13})(:,clip:end-clip),...
%     fld_names{15}, radar_tmp.(fld_names{15})(:,clip:end-clip,:));
% 
% radar_25km.(fld_names{16}) = radar_tmp.(fld_names{16})(clip:end-clip);
% radar_25km.(fld_names{17}) = radar_tmp.(fld_names{17})(clip:end-clip);
% 
% % Save file in output directory
% flnm = 'OIB_25km.mat';
% output_path = fullfile(data_path, output_dir, flnm);
% save(output_path, '-struct', 'radar_25km', '-v7.3')
% 
% % Clear saved variables (reduces memory load on additional datasets)
% clearvars -except data_path OIB_dir files output_dir cores Ndraw
% toc
% Processing times: [11.6s 267s 284s] 

radar_25km = load(fullfile(data_path, ...
    'radar/SEAT_Traverses/results_data/comparisons/dataset_length/OIB_25km.mat'));

%% radar_50km

% tic
% % Initialize for loop with first file in directory
% data_struct = OIB_import(fullfile(OIB_dir, files(7).name));
% fld_names = fieldnames(data_struct);
% data_cells = struct2cell(data_struct);
% 
% for i = 8:16
%     data_i = struct2cell(OIB_import(fullfile(OIB_dir, files(i).name)));
%     data_cells = cellfun(@horzcat, data_cells, data_i, 'UniformOutput', 0);
% end
% 
% % Average values for TWTT, time_trace, and collect_date
% data_cells{8} = median(data_cells{8}, 2);
% data_cells{9} = median(data_cells{9});
% data_cells{10} = median(data_cells{10});
% 
% % Flip data (so that the processing direction matches SEAT direction from
% % SEAT10-4 to SEAT10-6)
% data_flip = cellfun(@fliplr, data_cells, 'UniformOutput', false);
% 
% radar_tmp = cell2struct(data_flip, fld_names, 1);
% radar_tmp.dist = pathdist(radar_tmp.lat, radar_tmp.lon);
% toc
% 
% [radar_tmp] = radar_RT(radar_tmp, cores, Ndraw);
% [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
% toc
% 
% clip = 100;
% fld_names = fieldnames(radar_tmp);
% radar_50km = struct(fld_names{1}, radar_tmp.(fld_names{1}), fld_names{2}, ...
%     radar_tmp.(fld_names{2})(clip:end-clip), ...
%     fld_names{3}, radar_tmp.(fld_names{3})(clip:end-clip),...
%     fld_names{4}, radar_tmp.(fld_names{4})(clip:end-clip),...
%     fld_names{5}, radar_tmp.(fld_names{5}), ...
%     fld_names{9}, radar_tmp.(fld_names{9})(:,clip:end-clip),...
%     fld_names{10}, radar_tmp.(fld_names{10})(:,clip:end-clip),...
%     fld_names{11}, radar_tmp.(fld_names{11})(:,clip:end-clip),...
%     fld_names{13}, radar_tmp.(fld_names{13})(:,clip:end-clip),...
%     fld_names{15}, radar_tmp.(fld_names{15})(:,clip:end-clip,:));
% 
% radar_50km.(fld_names{16}) = radar_tmp.(fld_names{16})(clip:end-clip);
% radar_50km.(fld_names{17}) = radar_tmp.(fld_names{17})(clip:end-clip);
% 
% % Save file in output directory
% flnm = 'OIB_50km.mat';
% output_path = fullfile(data_path, output_dir, flnm);
% save(output_path, '-struct', 'radar_50km', '-v7.3')
% 
% % Clear saved variables (reduces memory load on additional datasets)
% clearvars -except data_path OIB_dir files output_dir cores Ndraw
% toc
% Processing times: [24s 764s 783s] 

radar_50km = load(fullfile(data_path, ...
    'radar/SEAT_Traverses/results_data/comparisons/dataset_length/OIB_50km.mat'));

%% radar_100km

% tic
% % Initialize for loop with first file in directory
% data_struct = OIB_import(fullfile(OIB_dir, files(1).name));
% fld_names = fieldnames(data_struct);
% data_cells = struct2cell(data_struct);
% 
% for i = 2:20
%     data_i = struct2cell(OIB_import(fullfile(OIB_dir, files(i).name)));
%     data_cells = cellfun(@horzcat, data_cells, data_i, 'UniformOutput', 0);
% end
% 
% % Average values for TWTT, time_trace, and collect_date
% data_cells{8} = median(data_cells{8}, 2);
% data_cells{9} = median(data_cells{9});
% data_cells{10} = median(data_cells{10});
% 
% % Flip data (so that the processing direction matches SEAT direction from
% % SEAT10-4 to SEAT10-6)
% data_flip = cellfun(@fliplr, data_cells, 'UniformOutput', false);
% 
% radar_tmp = cell2struct(data_flip, fld_names, 1);
% radar_tmp.dist = pathdist(radar_tmp.lat, radar_tmp.lon);
% toc
% 
% [radar_tmp] = radar_RT(radar_tmp, cores, Ndraw);
% [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
% toc
% 
% clip = 100;
% fld_names = fieldnames(radar_tmp);
% radar_100km = struct(fld_names{1}, radar_tmp.(fld_names{1}), fld_names{2}, ...
%     radar_tmp.(fld_names{2})(clip:end-clip), ...
%     fld_names{3}, radar_tmp.(fld_names{3})(clip:end-clip),...
%     fld_names{4}, radar_tmp.(fld_names{4})(clip:end-clip),...
%     fld_names{5}, radar_tmp.(fld_names{5}), ...
%     fld_names{9}, radar_tmp.(fld_names{9})(:,clip:end-clip),...
%     fld_names{10}, radar_tmp.(fld_names{10})(:,clip:end-clip),...
%     fld_names{11}, radar_tmp.(fld_names{11})(:,clip:end-clip),...
%     fld_names{13}, radar_tmp.(fld_names{13})(:,clip:end-clip),...
%     fld_names{15}, radar_tmp.(fld_names{15})(:,clip:end-clip,:));
% 
% radar_100km.(fld_names{16}) = radar_tmp.(fld_names{16})(clip:end-clip);
% radar_100km.(fld_names{17}) = radar_tmp.(fld_names{17})(clip:end-clip);
% 
% % Save file in output directory
% flnm = 'OIB_100km.mat';
% output_path = fullfile(data_path, output_dir, flnm);
% save(output_path, '-struct', 'radar_100km', '-v7.3')
% 
% % Clear saved variables (reduces memory load on additional datasets)
% clearvars -except data_path OIB_dir files output_dir cores Ndraw
% toc
% Processing times: [49s 2450s 2509s] 

radar_100km = load(fullfile(data_path, ...
    'radar/SEAT_Traverses/results_data/comparisons/dataset_length/OIB_100km.mat'));

%%

radar_full = radar_100km;

bias = [];

dist_50km = pdist2([radar_50km.Easting' radar_50km.Northing'], ...
    [radar_full.Easting' radar_full.Northing']);
[~, dist_idx] = min(dist_50km, [], 2);
ages_sub = median(radar_50km.ages, 3);
ages_full = median(radar_full.ages(:,dist_idx,:), 3);
bias = [bias median(ages_full(end,:) - ages_sub(end,:))];

figure(1)
hold on
h1 = histogram(ages_full(end,:)-ages_sub(end,:), 50);

figure(2)
hold on
plot(ages_full(end,:) - ages_sub(end,:), '.')

dist_25km = pdist2([radar_25km.Easting' radar_25km.Northing'], ...
    [radar_full.Easting' radar_full.Northing']);
[~, dist_idx] = min(dist_25km, [], 2);
ages_sub = median(radar_25km.ages, 3);
ages_full = median(radar_full.ages(:,dist_idx,:), 3);
bias = [bias median(ages_full(end,:)-ages_sub(end,:))];

figure(1)
hold on
h2 = histogram(ages_full(end,:)-ages_sub(end,:), 50);

figure(2)
hold on
plot(ages_full(end,:) - ages_sub(end,:), '.')

dist_15km = pdist2([radar_15km.Easting' radar_15km.Northing'], ...
    [radar_full.Easting' radar_full.Northing']);
[~, dist_idx] = min(dist_15km, [], 2);
ages_sub = median(radar_15km.ages, 3);
ages_full = median(radar_full.ages(:,dist_idx,:), 3);
bias = [bias median(ages_full(end,:)-ages_sub(end,:))];

figure(1)
hold on
h3 = histogram(ages_full(end,:)-ages_sub(end,:), 50);
hold off

figure(2)
hold on
plot(ages_full(end,:) - ages_sub(end,:), '.')
hold off
