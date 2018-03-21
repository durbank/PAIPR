% Script for comparing accumulation results from various different records
% (firn cores, SEAT Ku radar, OIB Ku radar, OIB accum radar, etc.)

%% Wrapper head

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/WAIS Variability/';
        addon_path = 'E:/Research/Antarctica/WAIS Variability/Addons/';
    case false
        data_path = '/Volumes/WARP/Research/Antarctica/WAIS Variability/';
        addon_path = '/Users/Durbank/Documents/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Import firn core data (eventually will incorporate this into a .mat file
% instead of excel file)
[cores] = import_cores(strcat(data_path, ['SEAT_cores' filesep ...
    'DGK_core_data.xlsx']));

%% Define radar files to import/process

% SEAT10-1
OIB_SNO_1_2011 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_211.nc';
OIB_SNO_1_2016 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2016/IRSNO1B_20161109_02_320.nc';
OIB_KU_1_2016 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Kuband/2016/IRKUB1B_20161109_02_320.nc';
SEAT10_1_files = {OIB_SNO_1_2011 OIB_SNO_1_2016 OIB_KU_1_2016};

% SEAT10-4
SEAT_KU_4_2010 = '/Volumes/WARP/Research/Antarctica/WAIS Variability/SEAT_Traverses/SEAT2010Kuband/ProcessedSEAT2010/grid_SEAT10_4/layers_ku_band_12132010_out_dec02.mat';
OIB_SNO_4_2011 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_272.nc';
OIB_SNO_4_2016 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2016/IRSNO1B_20161109_02_381.nc';
OIB_KU_4_2016 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Kuband/2016/IRKUB1B_20161109_02_381.nc';
SEAT10_4_files = {SEAT_KU_4_2010 OIB_SNO_4_2011 OIB_SNO_4_2016 OIB_KU_4_2016};

% SEAT10-5
SEAT_KU_5_2010 = '/Volumes/WARP/Research/Antarctica/WAIS Variability/SEAT_Traverses/SEAT2010Kuband/ProcessedSEAT2010/grid_SEAT10_5/layers_Ku_band_12152010_out_dec04.mat';
OIB_SNO_5_2011 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_257.nc';
OIB_SNO_5_2016 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2016/IRSNO1B_20161109_02_366.nc';
OIB_KU_5_2016 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Kuband/2016/IRKUB1B_20161109_02_366.nc';
SEAT10_5_files = {SEAT_KU_5_2010 OIB_SNO_5_2011 OIB_SNO_5_2016 OIB_KU_5_2016};

% SEAT10-6
SEAT_KU_6_2010 = '/Volumes/WARP/Research/Antarctica/WAIS Variability/SEAT_Traverses/SEAT2010Kuband/ProcessedSEAT2010/grid_SEAT10_6/layers_Ku_band_12162010_out_dec02.mat';
OIB_SNOW_6_2011 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_242.nc';
OIB_SNO_6_2016 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2016/IRSNO1B_20161109_02_350.nc';
OIB_KU_6_2016 = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Kuband/2016/IRKUB1B_20161109_02_350.nc';
SEAT10_6_files = {SEAT_KU_6_2010 OIB_SNOW_6_2011 OIB_SNO_6_2016 OIB_KU_6_2016};


%% Process radar data

Ndraw = 100;

%%SEAT10_1

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_1_files{1}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_1.OIB2011_SNO, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_1_files{2}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_1.OIB2016_SNO, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_1_files{3}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_1.OIB2016_KU, ~] = calc_SWE(radar, core, Ndraw);

%%SEAT10_4

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_4_files{1}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_4.SEAT_KU, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_4_files{2}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_4.OIB2011_SNO, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_4_files{3}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_4.OIB2016_SNO, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_4_files{4}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_4.OIB2016_KU, ~] = calc_SWE(radar, core, Ndraw);

%%SEAT10_5

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_5_files{1}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_5.SEAT_KU, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_5_files{2}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_5.OIB2011_SNO, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_5_files{3}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_5.OIB2016_SNO, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_5_files{4}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_5.OIB2016_KU, ~] = calc_SWE(radar, core, Ndraw);

%%SEAT10_6

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_6_files{1}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_6.SEAT_KU, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_6_files{2}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_6.OIB2011_SNO, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_6_files{3}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_6.OIB2016_SNO, ~] = calc_SWE(radar, core, Ndraw);

% Calculate radar ages and associated other data
[radar, core] = radar_age(SEAT10_6_files{4}, cores, Ndraw);
% Calculate annual accumulation rates from data
[SEAT10_6.OIB2016_KU, ~] = calc_SWE(radar, core, Ndraw);

%% Calculate core accumulations for each site

% Calculate accumulation for each core
core_name = {'SEAT10_1' 'SEAT10_4', 'SEAT10_5', 'SEAT10_6'};

for i = core_name
    core = cores.(i{1});
    core_accum_dt = 0.02*(1000*core.rho);
    yr_top = floor(core.age(2));
    yr_end = ceil(core.age(end));
    core_yr = (yr_top:-1:yr_end)';
    core_yr_idx = logical([1; diff(floor(core.age))]);
    yr_loc = find(core_yr_idx);
    
    core_accum = zeros(length(core_yr), Ndraw);
    for j = 1:Ndraw
        yr_loc_j = yr_loc;
        yr_loc_j(2:end-1) = yr_loc(2:end-1) + ...
            round(1*(mean(diff(yr_loc))/12)*randn(length(yr_loc)-2, 1));
        loc_idx = yr_loc_j<1;
        yr_loc_j(loc_idx) = yr_loc(loc_idx);
        
        core_accum_j = zeros(length(core_yr), 1);
        for n = 1:length(core_yr)
            core_accum_j(n) = sum(core_accum_dt(yr_loc_j(n)+1:yr_loc_j(n+1)));
        end
        core_accum(:,j) = core_accum_j;
    end
    cores.(i{1}).SMB_yr = core_yr;
    cores.(i{1}).SMB = core_accum;
    
end

% Add cores to SEAT structures
SEAT10_1.core = cores.SEAT10_1;
SEAT10_4.core = cores.SEAT10_4;
SEAT10_5.core = cores.SEAT10_5;
SEAT10_6.core = cores.SEAT10_6;

%% Find traces nearest to the relevent core

%%SEAT10-1
fdn = fieldnames(SEAT10_1);
for i = 1:length(fdn)-1
    
    [distance] = pdist2([SEAT10_1.core.Easting SEAT10_1.core.Northing],...
        [SEAT10_1.(fdn{i}).Easting' SEAT10_1.(fdn{i}).Northing'], 'Euclidean');
    [~, idx_near.SEAT10_1.(fdn{i})] = min(distance);
end

%%SEAT10-4
fdn = fieldnames(SEAT10_4);
for i = 1:length(fdn)-1
    
    [distance] = pdist2([SEAT10_4.core.Easting SEAT10_4.core.Northing],...
        [SEAT10_4.(fdn{i}).Easting' SEAT10_4.(fdn{i}).Northing'], 'Euclidean');
    [~, idx_near.SEAT10_4.(fdn{i})] = min(distance);
end

%%SEAT10-5
fdn = fieldnames(SEAT10_5);
for i = 1:length(fdn)-1
    
    [distance] = pdist2([SEAT10_5.core.Easting SEAT10_5.core.Northing],...
        [SEAT10_5.(fdn{i}).Easting' SEAT10_5.(fdn{i}).Northing'], 'Euclidean');
    [~, idx_near.SEAT10_5.(fdn{i})] = min(distance);
end

%%SEAT10-6
fdn = fieldnames(SEAT10_6);
for i = 1:length(fdn)-1
    
    [distance] = pdist2([SEAT10_6.core.Easting SEAT10_6.core.Northing],...
        [SEAT10_6.(fdn{i}).Easting' SEAT10_6.(fdn{i}).Northing'], 'Euclidean');
    [~, idx_near.SEAT10_6.(fdn{i})] = min(distance);
end



%% Diagnostic figures 

% load('/Users/Durbank/Documents/MATLAB/Research/tmp_data/NASA_plots.mat')

%%SEAT10_4 plots

% Plot SEAT Ku radargram
figure('Position', [200 200 1000 500])
imagesc(SEAT10_4.SEAT_KU.dist, SEAT10_4.SEAT_KU.depth, ...
    SEAT10_4.SEAT_KU.data_smooth, [-2 2])
colorbar
title('SEAT Ku radargram for Site SEAT10-4')
xlabel('Distance along profile (m)')
ylabel('Depth (m)')

% Plot OIB snow radargram
figure('Position', [200 200 1000 500])
imagesc(SEAT10_4.OIB2011_SNO.dist, SEAT10_4.OIB2011_SNO.depth, ...
    SEAT10_4.OIB2011_SNO.data_smooth, [-2 2])
colorbar
title('OIB-snow radargram for Site SEAT10-4')
xlabel('Distance along profile (m)')
ylabel('Depth (m)')

% Plot age-depth profiles
figure('Position', [200 200 1000 500])
hold on
h1 = plot(SEAT10_4.core.depth, SEAT10_4.core.age, 'b', 'LineWidth', 2);
h2 = plot(SEAT10_4.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_4.SEAT_KU.age(:,idx_near.SEAT10_4.SEAT_KU,:)), 2),...
    'r', 'LineWidth', 2);
plot(SEAT10_4.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_4.SEAT_KU.age(:,idx_near.SEAT10_4.SEAT_KU,:)), 2)...
    + 2*std(squeeze(SEAT10_4.SEAT_KU.age(:,idx_near.SEAT10_4.SEAT_KU,:)),[],2),...
    'r--');
plot(SEAT10_4.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_4.SEAT_KU.age(:,idx_near.SEAT10_4.SEAT_KU,:)), 2)...
    - 2*std(squeeze(SEAT10_4.SEAT_KU.age(:,idx_near.SEAT10_4.SEAT_KU,:)),[],2),...
    'r--');
h3 = plot(SEAT10_4.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_4.OIB2011_SNO.age(:,idx_near.SEAT10_4.OIB2011_SNO,:)), 2),...
    'c', 'LineWidth', 2);
plot(SEAT10_4.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_4.OIB2011_SNO.age(:,idx_near.SEAT10_4.OIB2011_SNO,:)), 2)...
    + 2*std(squeeze(SEAT10_4.OIB2011_SNO.age(:,idx_near.SEAT10_4.OIB2011_SNO,:)),[],2),...
    'c--');
plot(SEAT10_4.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_4.OIB2011_SNO.age(:,idx_near.SEAT10_4.OIB2011_SNO,:)), 2)...
    - 2*std(squeeze(SEAT10_4.OIB2011_SNO.age(:,idx_near.SEAT10_4.OIB2011_SNO,:)),[],2),...
    'c--');
legend([h1 h2 h3], 'Firn core', 'SEAT Ku radar', 'OIB snow radar')
xlabel('Depth (meters)')
ylabel('Calendar year')
title('Comparison of site SEAT10-4 age-depth profiles')
hold off

% Plot annual accumulation
figure('Position', [200 200 1000 500])
hold on
h1 = plot(SEAT10_4.core.SMB_yr, mean(SEAT10_4.core.SMB, 2), 'b', 'LineWidth', 2);
plot(SEAT10_4.core.SMB_yr, ...
    mean(SEAT10_4.core.SMB, 2)+2*std(SEAT10_4.core.SMB, [], 2), 'b--')
plot(SEAT10_4.core.SMB_yr, ...
    mean(SEAT10_4.core.SMB, 2)-2*std(SEAT10_4.core.SMB, [], 2), 'b--')

h2 = plot(SEAT10_4.SEAT_KU.SMB_yr{idx_near.SEAT10_4.SEAT_KU}, ...
    mean(SEAT10_4.SEAT_KU.SMB{idx_near.SEAT10_4.SEAT_KU}, 2), 'r', 'LineWidth', 2);
plot(SEAT10_4.SEAT_KU.SMB_yr{idx_near.SEAT10_4.SEAT_KU}, ...
    mean(SEAT10_4.SEAT_KU.SMB{idx_near.SEAT10_4.SEAT_KU}, 2) + ...
    std(SEAT10_4.SEAT_KU.SMB{idx_near.SEAT10_4.SEAT_KU}, [], 2), 'r--')
plot(SEAT10_4.SEAT_KU.SMB_yr{idx_near.SEAT10_4.SEAT_KU}, ...
    mean(SEAT10_4.SEAT_KU.SMB{idx_near.SEAT10_4.SEAT_KU}, 2) - ...
    std(SEAT10_4.SEAT_KU.SMB{idx_near.SEAT10_4.SEAT_KU}, [], 2), 'r--')

h3 = plot(SEAT10_4.OIB2011_SNO.SMB_yr{idx_near.SEAT10_4.OIB2011_SNO}, ...
    mean(SEAT10_4.OIB2011_SNO.SMB{idx_near.SEAT10_4.OIB2011_SNO}, 2), 'c', 'LineWidth', 2);
plot(SEAT10_4.OIB2011_SNO.SMB_yr{idx_near.SEAT10_4.OIB2011_SNO}, ...
    mean(SEAT10_4.OIB2011_SNO.SMB{idx_near.SEAT10_4.OIB2011_SNO}, 2) + ...
    2*std(SEAT10_4.OIB2011_SNO.SMB{idx_near.SEAT10_4.OIB2011_SNO}, [], 2), 'c--')
plot(SEAT10_4.OIB2011_SNO.SMB_yr{idx_near.SEAT10_4.OIB2011_SNO}, ...
    mean(SEAT10_4.OIB2011_SNO.SMB{idx_near.SEAT10_4.OIB2011_SNO}, 2) - ...
    2*std(SEAT10_4.OIB2011_SNO.SMB{idx_near.SEAT10_4.OIB2011_SNO}, [], 2), 'c--')

legend([h1 h2 h3], 'Firn core', 'SEAT Ku radar', 'OIB snow radar')
xlabel('Calendar Year')
ylabel('Accumulation (mm w.e.)')
title('Comparison of site SEAT10-4 annual accumulation')
hold off

%%SEAT10_5 plots

% Plot SEAT Ku radargram
figure('Position', [200 200 1000 500])
imagesc(SEAT10_5.SEAT_KU.dist, SEAT10_5.SEAT_KU.depth, ...
    SEAT10_5.SEAT_KU.data_smooth, [-2 2])
colorbar
title('SEAT Ku radargram for Site SEAT10-5')
xlabel('Distance along profile (m)')
ylabel('Depth (m)')

% Plot OIB snow radargram
figure('Position', [200 200 1000 500])
imagesc(SEAT10_5.OIB2011_SNO.dist, SEAT10_5.OIB2011_SNO.depth, ...
    SEAT10_5.OIB2011_SNO.data_smooth, [-2 2])
colorbar
title('OIB-snow radargram for Site SEAT10-5')
xlabel('Distance along profile (m)')
ylabel('Depth (m)')

% Plot age-depth profiles
figure('Position', [200 200 1000 500])
hold on
h1 = plot(SEAT10_5.core.depth, SEAT10_5.core.age, 'b', 'LineWidth', 2);
h2 = plot(SEAT10_5.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_5.SEAT_KU.age(:,idx_near.SEAT10_5.SEAT_KU,:)), 2),...
    'r', 'LineWidth', 2);
plot(SEAT10_5.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_5.SEAT_KU.age(:,idx_near.SEAT10_5.SEAT_KU,:)), 2)...
    + 2*std(squeeze(SEAT10_5.SEAT_KU.age(:,idx_near.SEAT10_5.SEAT_KU,:)),[],2),...
    'r--');
plot(SEAT10_5.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_5.SEAT_KU.age(:,idx_near.SEAT10_5.SEAT_KU,:)), 2)...
    - 2*std(squeeze(SEAT10_5.SEAT_KU.age(:,idx_near.SEAT10_5.SEAT_KU,:)),[],2),...
    'r--');
h3 = plot(SEAT10_5.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_5.OIB2011_SNO.age(:,idx_near.SEAT10_5.OIB2011_SNO,:)), 2),...
    'c', 'LineWidth', 2);
plot(SEAT10_5.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_5.OIB2011_SNO.age(:,idx_near.SEAT10_5.OIB2011_SNO,:)), 2)...
    + 2*std(squeeze(SEAT10_5.OIB2011_SNO.age(:,idx_near.SEAT10_5.OIB2011_SNO,:)),[],2),...
    'c--');
plot(SEAT10_5.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_5.OIB2011_SNO.age(:,idx_near.SEAT10_5.OIB2011_SNO,:)), 2)...
    - 2*std(squeeze(SEAT10_5.OIB2011_SNO.age(:,idx_near.SEAT10_5.OIB2011_SNO,:)),[],2),...
    'c--');
legend([h1 h2 h3], 'Firn core', 'SEAT Ku radar', 'OIB snow radar')
xlabel('Depth (meters)')
ylabel('Calendar year')
title('Comparison of site SEAT10-5 age-depth profiles')
hold off

% Plot annual accumulation
figure('Position', [200 200 1000 500])
hold on
h1 = plot(SEAT10_5.core.SMB_yr, mean(SEAT10_5.core.SMB, 2), 'b', 'LineWidth', 2);
plot(SEAT10_5.core.SMB_yr, ...
    mean(SEAT10_5.core.SMB, 2) + 2*std(SEAT10_5.core.SMB, [], 2), 'b--')
plot(SEAT10_5.core.SMB_yr, ...
    mean(SEAT10_5.core.SMB, 2) - 2*std(SEAT10_5.core.SMB, [], 2), 'b--')

h2 = plot(SEAT10_5.SEAT_KU.SMB_yr{idx_near.SEAT10_5.SEAT_KU}, ...
    mean(SEAT10_5.SEAT_KU.SMB{idx_near.SEAT10_5.SEAT_KU}, 2), 'r', 'LineWidth', 2);
plot(SEAT10_5.SEAT_KU.SMB_yr{idx_near.SEAT10_5.SEAT_KU}, ...
    mean(SEAT10_5.SEAT_KU.SMB{idx_near.SEAT10_5.SEAT_KU}, 2) + ...
    2*std(SEAT10_5.SEAT_KU.SMB{idx_near.SEAT10_5.SEAT_KU}, [], 2), 'r--')
plot(SEAT10_5.SEAT_KU.SMB_yr{idx_near.SEAT10_5.SEAT_KU}, ...
    mean(SEAT10_5.SEAT_KU.SMB{idx_near.SEAT10_5.SEAT_KU}, 2) - ...
    2*std(SEAT10_5.SEAT_KU.SMB{idx_near.SEAT10_5.SEAT_KU}, [], 2), 'r--')

h3 = plot(SEAT10_5.OIB2011_SNO.SMB_yr{idx_near.SEAT10_5.OIB2011_SNO}, ...
    mean(SEAT10_5.OIB2011_SNO.SMB{idx_near.SEAT10_5.OIB2011_SNO}, 2), 'c', 'LineWidth', 2);
plot(SEAT10_5.OIB2011_SNO.SMB_yr{idx_near.SEAT10_5.OIB2011_SNO}, ...
    mean(SEAT10_5.OIB2011_SNO.SMB{idx_near.SEAT10_5.OIB2011_SNO}, 2) + ...
    2*std(SEAT10_5.OIB2011_SNO.SMB{idx_near.SEAT10_5.OIB2011_SNO}, [], 2), 'c--')
plot(SEAT10_5.OIB2011_SNO.SMB_yr{idx_near.SEAT10_5.OIB2011_SNO}, ...
    mean(SEAT10_5.OIB2011_SNO.SMB{idx_near.SEAT10_5.OIB2011_SNO}, 2) - ...
    2*std(SEAT10_5.OIB2011_SNO.SMB{idx_near.SEAT10_5.OIB2011_SNO}, [], 2), 'c--')

legend([h1 h2 h3], 'Firn core', 'SEAT Ku radar', 'OIB snow radar')
xlabel('Calendar Year')
ylabel('Accumulation (mm w.e.)')
title('Comparison of site SEAT10-5 annual accumulation')
hold off

%%SEAT10_6 plots

% Plot SEAT Ku radargram
figure('Position', [200 200 1000 500])
imagesc(SEAT10_6.SEAT_KU.dist, SEAT10_6.SEAT_KU.depth, ...
    SEAT10_6.SEAT_KU.data_smooth, [-2 2])
colorbar
title('SEAT Ku radargram for Site SEAT10-6')
xlabel('Distance along profile (m)')
ylabel('Depth (m)')

% Plot OIB snow radargram
figure('Position', [200 200 1000 500])
imagesc(SEAT10_6.OIB2011_SNO.dist, SEAT10_6.OIB2011_SNO.depth, ...
    SEAT10_6.OIB2011_SNO.data_smooth, [-2 2])
colorbar
title('OIB-snow radargram for Site SEAT10-6')
xlabel('Distance along profile (m)')
ylabel('Depth (m)')

% Plot age-depth profiles
figure('Position', [200 200 1000 500])
hold on
h1 = plot(SEAT10_6.core.depth, SEAT10_6.core.age, 'b', 'LineWidth', 2);
h2 = plot(SEAT10_6.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_6.SEAT_KU.age(:,idx_near.SEAT10_6.SEAT_KU,:)), 2),...
    'r', 'LineWidth', 2);
plot(SEAT10_6.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_6.SEAT_KU.age(:,idx_near.SEAT10_6.SEAT_KU,:)), 2)...
    + 2*std(squeeze(SEAT10_6.SEAT_KU.age(:,idx_near.SEAT10_6.SEAT_KU,:)),[],2),...
    'r--');
plot(SEAT10_6.SEAT_KU.depth, ...
    mean(squeeze(SEAT10_6.SEAT_KU.age(:,idx_near.SEAT10_6.SEAT_KU,:)), 2)...
    - 2*std(squeeze(SEAT10_6.SEAT_KU.age(:,idx_near.SEAT10_6.SEAT_KU,:)),[],2),...
    'r--');
h3 = plot(SEAT10_6.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_6.OIB2011_SNO.age(:,idx_near.SEAT10_6.OIB2011_SNO,:)), 2),...
    'c', 'LineWidth', 2);
plot(SEAT10_6.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_6.OIB2011_SNO.age(:,idx_near.SEAT10_6.OIB2011_SNO,:)), 2)...
    + 2*std(squeeze(SEAT10_6.OIB2011_SNO.age(:,idx_near.SEAT10_6.OIB2011_SNO,:)),[],2),...
    'c--');
plot(SEAT10_6.OIB2011_SNO.depth, ...
    mean(squeeze(SEAT10_6.OIB2011_SNO.age(:,idx_near.SEAT10_6.OIB2011_SNO,:)), 2)...
    - 2*std(squeeze(SEAT10_6.OIB2011_SNO.age(:,idx_near.SEAT10_6.OIB2011_SNO,:)),[],2),...
    'c--');
legend([h1 h2 h3], 'Firn core', 'SEAT Ku radar', 'OIB snow radar')
xlabel('Depth (meters)')
ylabel('Calendar year')
title('Comparison of site SEAT10-6 age-depth profiles')
hold off

% Plot annual accumulation
figure('Position', [200 200 1000 500])
hold on
h1 = plot(SEAT10_6.core.SMB_yr, mean(SEAT10_6.core.SMB, 2), 'b', 'LineWidth', 2);
plot(SEAT10_6.core.SMB_yr, ...
    mean(SEAT10_6.core.SMB, 2) + 2*std(SEAT10_6.core.SMB, [], 2), 'b--')
plot(SEAT10_6.core.SMB_yr, ...
    mean(SEAT10_6.core.SMB, 2) - 2*std(SEAT10_6.core.SMB, [], 2), 'b--')

h2 = plot(SEAT10_6.SEAT_KU.SMB_yr{idx_near.SEAT10_6.SEAT_KU}, ...
    mean(SEAT10_6.SEAT_KU.SMB{idx_near.SEAT10_6.SEAT_KU}, 2), 'r', 'LineWidth', 2);
plot(SEAT10_6.SEAT_KU.SMB_yr{idx_near.SEAT10_6.SEAT_KU}, ...
    mean(SEAT10_6.SEAT_KU.SMB{idx_near.SEAT10_6.SEAT_KU}, 2) + ...
    2*std(SEAT10_6.SEAT_KU.SMB{idx_near.SEAT10_6.SEAT_KU}, [], 2), 'r--')
plot(SEAT10_6.SEAT_KU.SMB_yr{idx_near.SEAT10_6.SEAT_KU}, ...
    mean(SEAT10_6.SEAT_KU.SMB{idx_near.SEAT10_6.SEAT_KU}, 2) - ...
    2*std(SEAT10_6.SEAT_KU.SMB{idx_near.SEAT10_6.SEAT_KU}, [], 2), 'r--')

h3 = plot(SEAT10_6.OIB2011_SNO.SMB_yr{idx_near.SEAT10_6.OIB2011_SNO}, ...
    mean(SEAT10_6.OIB2011_SNO.SMB{idx_near.SEAT10_6.OIB2011_SNO}, 2), 'c', 'LineWidth', 2);
plot(SEAT10_6.OIB2011_SNO.SMB_yr{idx_near.SEAT10_6.OIB2011_SNO}, ...
    mean(SEAT10_6.OIB2011_SNO.SMB{idx_near.SEAT10_6.OIB2011_SNO}, 2) + ...
    2*std(SEAT10_6.OIB2011_SNO.SMB{idx_near.SEAT10_6.OIB2011_SNO}, [], 2), 'c--')
plot(SEAT10_6.OIB2011_SNO.SMB_yr{idx_near.SEAT10_6.OIB2011_SNO}, ...
    mean(SEAT10_6.OIB2011_SNO.SMB{idx_near.SEAT10_6.OIB2011_SNO}, 2) - ...
    2*std(SEAT10_6.OIB2011_SNO.SMB{idx_near.SEAT10_6.OIB2011_SNO}, [], 2), 'c--')

legend([h1 h2 h3], 'Firn core', 'SEAT Ku radar', 'OIB snow radar')
xlabel('Calendar Year')
ylabel('Accumulation (mm w.e.)')
title('Comparison of site SEAT10-6 annual accumulation')
hold off

clearvars -except cores idx_near SEAT10_1 SEAT10_4 SEAT10_5 SEAT10_6
