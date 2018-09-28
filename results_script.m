% Script for assessing and comparing results of accum-radar, and generating
% figures for the methods paper

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
% Add CReSIS OIB MATLAB reader functions to path
addon_folder = fullfile(addon_path, 'cresis-L1B-matlab-readers/');
addpath(genpath(addon_folder))

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Define number of Monte Carlo simulations to perform
Ndraw = 100;

%%

wild = '*.mat';
SEAT_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
    'SEAT2010Kuband/SEAT10_4toSEAT10_6/SMB_results/', wild));

SEAT_E = [];
SEAT_N = [];
SEAT_SMB_MC = [];
SEAT_yr = [];
for i = 1:length(SEAT_files)
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'Easting');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'Northing');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'SMB');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'SMB_yr');
    
    SEAT_E = [SEAT_E Easting];
    SEAT_N = [SEAT_N Northing];
    SEAT_SMB_MC = [SEAT_SMB_MC SMB];
    SEAT_yr = [SEAT_yr SMB_yr];
end
clear Easting Northing SMB SMB_yr

SEAT_SMB = cellfun(@(x) median(x, 2), SEAT_SMB_MC, 'UniformOutput', 0);
SEAT_std = cellfun(@(x) std(x, [], 2), SEAT_SMB_MC, 'UniformOutput', 0);


wild = '*.mat';
OIB_files = dir(fullfile(data_path, ...
    'IceBridge/SEAT10_4to10_6/2011_SNO/SMB_results', wild));

OIB_E = [];
OIB_N = [];
OIB_SMB = [];
OIB_yr = [];
for i = 1:length(OIB_files)
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Easting');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Northing');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB_yr');
    
    OIB_E = [OIB_E Easting];
    OIB_N = [OIB_N Northing];
    OIB_SMB = [OIB_SMB SMB];
    OIB_yr = [OIB_yr SMB_yr];
end
clear Easting Northing SMB SMB_yr i wild


near_dist = 10000;


%% Compare SMB estimates within 10 km of site 4 (core site SEAT10-4)


D_SEAT4 = pdist2([cores.SEAT10_4.Easting, cores.SEAT10_4.Northing],...
    [SEAT_E' SEAT_N']);
SEAT4_idx = D_SEAT4 <= near_dist;

D_OIB4 = pdist2([cores.SEAT10_4.Easting, cores.SEAT10_4.Northing],...
    [OIB_E' OIB_N']);
OIB4_idx = D_OIB4 <= near_dist;

yr_start = min([cores.SEAT10_4.SMB_yr(1) cellfun(@(x) x(1), ...
    SEAT_yr(SEAT4_idx)) cellfun(@(x) x(1), OIB_yr(OIB4_idx))]);
yr_end = max([cores.SEAT10_4.SMB_yr(end) cellfun(@(x) x(end), ...
    SEAT_yr(SEAT4_idx)) cellfun(@(x) x(end), OIB_yr(OIB4_idx))]);
site4_yr = (yr_start:-1:yr_end)';

SEAT_start = cellfun(@(x) find(x==yr_start, 1), SEAT_yr(SEAT4_idx));
SEAT_end = cellfun(@(x) find(x==yr_end, 1), SEAT_yr(SEAT4_idx));
OIB_start = cellfun(@(x) find(x==yr_start, 1), OIB_yr(OIB4_idx));
OIB_end = cellfun(@(x) find(x==yr_end, 1), OIB_yr(OIB4_idx));

SMB_tmp = SEAT_SMB(SEAT4_idx);
SEAT4_SMB = zeros(length(site4_yr), length(SEAT_start));
for i = 1:length(SEAT_start)
    
    SEAT4_SMB(:,i) = SMB_tmp{i}(SEAT_start(i):SEAT_end(i));
end
SEAT4_SMB = movmean(SEAT4_SMB, 3);

SMB_tmp = OIB_SMB(OIB4_idx);
OIB4_SMB = zeros(length(site4_yr), length(OIB_start));
for i = 1:length(OIB_start)
    
    OIB4_SMB(:,i) = SMB_tmp{i}(OIB_start(i):OIB_end(i));
end
OIB4_SMB = movmean(OIB4_SMB, 3);

figure
hold on
for n = 1:size(SEAT4_SMB, 2)
    h0 = plot(site4_yr, SEAT4_SMB(:,n), 'r', 'LineWidth', 0.5);
    h0.Color(4) = 0.02;
end
for n = 1:size(OIB4_SMB, 2)
    h1 = plot(site4_yr, OIB4_SMB(:,n), 'm', 'LineWidth', 0.5);
    h1.Color(4) = 0.02;
end
core_start = find(cores.SEAT10_4.SMB_yr==yr_start);
core_end = find(cores.SEAT10_4.SMB_yr==yr_end);
core4_SMB = cores.SEAT10_4.SMB(core_start:core_end,:);
core4_SMB = movmean(core4_SMB, 3);
h2 = plot(cores.SEAT10_4.SMB_yr(core_start:core_end), median(core4_SMB, 2), ...
    'b', 'LineWidth', 2);
plot(cores.SEAT10_4.SMB_yr(core_start:core_end), median(core4_SMB, 2) + ...
    std(core4_SMB, [], 2), 'b--', 'LineWidth', 0.5);
plot(cores.SEAT10_4.SMB_yr(core_start:core_end), median(core4_SMB, 2) - ...
    std(core4_SMB, [], 2), 'b--', 'LineWidth', 0.5);
xlim([yr_end yr_start])
hold off



[p] = wanova(core4_SMB, SEAT4_SMB, OIB4_SMB);


%% Find SMB estimates within 10 km of SEAT10-5

D_SEAT5 = pdist2([cores.SEAT10_5.Easting, cores.SEAT10_5.Northing],...
    [SEAT_E' SEAT_N']);
SEAT5_idx = D_SEAT5 <= near_dist;

SEAT5_fullSMB = cellfun(@(x) median(x, 2), SEAT_SMB(SEAT5_idx), ...
    'UniformOutput', 0);
SEAT5_yr = SEAT_yr(SEAT5_idx);

D_OIB5 = pdist2([cores.SEAT10_5.Easting, cores.SEAT10_5.Northing],...
    [OIB_E' OIB_N']);
OIB5_idx = D_OIB5 <= near_dist;
OIB5_fullSMB = cellfun(@(x) median(x, 2), OIB_SMB(OIB5_idx), ...
    'UniformOutput', 0);
OIB5_yr = OIB_yr(OIB5_idx);


figure
hold on
for n = 1:length(SEAT5_yr)
    h0 = plot(SEAT5_yr{n}, SEAT5_fullSMB{n}, 'r', 'LineWidth', 0.5);
    h0.Color(4) = 0.02;
end
for n = 1:length(OIB5_yr)
    h1 = plot(OIB5_yr{n}, OIB5_fullSMB{n}, 'm', 'LineWidth', 0.5);
    h1.Color(4) = 0.02;
end
h2 = plot(cores.SEAT10_5.SMB_yr, median(cores.SEAT10_5.SMB, 2), ...
    'b', 'LineWidth', 2);
plot(cores.SEAT10_5.SMB_yr, median(cores.SEAT10_5.SMB, 2) + ...
    std(cores.SEAT10_5.SMB, [], 2), 'b--', 'LineWidth', 0.5);
plot(cores.SEAT10_5.SMB_yr, median(cores.SEAT10_5.SMB, 2) - ...
    std(cores.SEAT10_5.SMB, [], 2), 'b--', 'LineWidth', 0.5);
hold off

%% Find SMB estimates within 10 km of SEAT10-6

D_SEAT6 = pdist2([cores.SEAT10_6.Easting, cores.SEAT10_6.Northing],...
    [SEAT_E' SEAT_N']);
SEAT6_idx = D_SEAT6 <= near_dist;

SEAT6_fullSMB = cellfun(@(x) median(x, 2), SEAT_SMB(SEAT6_idx), ...
    'UniformOutput', 0);
SEAT6_yr = SEAT_yr(SEAT6_idx);

D_OIB6 = pdist2([cores.SEAT10_6.Easting, cores.SEAT10_6.Northing],...
    [OIB_E' OIB_N']);
OIB6_idx = D_OIB6 <= near_dist;
OIB6_fullSMB = cellfun(@(x) median(x, 2), OIB_SMB(OIB6_idx), ...
    'UniformOutput', 0);
OIB6_yr = OIB_yr(OIB6_idx);


figure
hold on
for n = 1:length(SEAT6_yr)
    h0 = plot(SEAT6_yr{n}, SEAT6_fullSMB{n}, 'r', 'LineWidth', 0.5);
    h0.Color(4) = 0.02;
end
for n = 1:length(OIB6_yr)
    h1 = plot(OIB6_yr{n}, OIB6_fullSMB{n}, 'm', 'LineWidth', 0.5);
    h1.Color(4) = 0.02;
end
h2 = plot(cores.SEAT10_6.SMB_yr, median(cores.SEAT10_6.SMB, 2), ...
    'b', 'LineWidth', 2);
plot(cores.SEAT10_6.SMB_yr, median(cores.SEAT10_6.SMB, 2) + ...
    std(cores.SEAT10_6.SMB, [], 2), 'b--', 'LineWidth', 0.5);
plot(cores.SEAT10_6.SMB_yr, median(cores.SEAT10_6.SMB, 2) - ...
    std(cores.SEAT10_6.SMB, [], 2), 'b--', 'LineWidth', 0.5);
hold off
