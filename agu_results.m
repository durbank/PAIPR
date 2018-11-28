% Script to perform analyses and generate figures for the results section
% of AGU2018 poster presentation

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'laptop';
        %         computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
                
            case 'laptop'
                data_path = 'F:/Research/Antarctica/Data/';
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

output_dir = uigetdir(data_path, ...
    'Select directory to which to output images');

%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);
Ndraw = 100;

% Find file names for previously processed SEAT2010 accumulation data
wild = '*.mat';
SEAT10_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
    'SEAT2010Kuband/RawSEAT2010/SMB_results/', wild));

% Preallocate arrays of sufficient size for data
seat10_E = zeros(1, length(SEAT10_files)*2*(50*1000/25));
seat10_N = seat10_E;
seat10_SMB_MC = cell(1, length(SEAT10_files)*2*(50*1000/25));
seat10_yr = seat10_SMB_MC;

for i = 1:length(SEAT10_files)
    
    % Load relevent data from current data file
    load(fullfile(SEAT10_files(i).folder, SEAT10_files(i).name), 'Easting');
    load(fullfile(SEAT10_files(i).folder, SEAT10_files(i).name), 'Northing');
    load(fullfile(SEAT10_files(i).folder, SEAT10_files(i).name), 'SMB');
    load(fullfile(SEAT10_files(i).folder, SEAT10_files(i).name), 'SMB_yr');
    
    % Find position of last data entered into preallocated arrays
    next_idx = sum(~cellfun(@isempty, seat10_SMB_MC)) + 1;
    
    % Fill current iteration data into preallocated arrays
    seat10_E(next_idx:next_idx+length(Easting)-1) = Easting;
    seat10_N(next_idx:next_idx+length(Northing)-1) = Northing;
    seat10_SMB_MC(next_idx:next_idx+length(SMB)-1) = SMB;
    seat10_yr(next_idx:next_idx+length(SMB_yr)-1) = SMB_yr;
    
end

% Find file names for previously processed SEAT2011 accumulation data
SEAT11_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
    'SEAT2011Kuband/RawSEAT2011/SMB_results/', wild));

% Preallocate arrays of sufficient size for data
seat11_E = zeros(1, length(SEAT11_files)*2*(50*1000/25));
seat11_N = seat11_E;
seat11_SMB_MC = cell(1, length(SEAT11_files)*2*(50*1000/25));
seat11_yr = seat10_SMB_MC;

for i = 1:length(SEAT11_files)
    
    % Load relevent data from current data file
    load(fullfile(SEAT11_files(i).folder, SEAT11_files(i).name), 'Easting');
    load(fullfile(SEAT11_files(i).folder, SEAT11_files(i).name), 'Northing');
    load(fullfile(SEAT11_files(i).folder, SEAT11_files(i).name), 'SMB');
    load(fullfile(SEAT11_files(i).folder, SEAT11_files(i).name), 'SMB_yr');
    
    % Find position of last data entered into preallocated arrays
    next_idx = sum(~cellfun(@isempty, seat11_SMB_MC)) + 1;
    
    % Fill current iteration data into preallocated arrays
    seat11_E(next_idx:next_idx+length(Easting)-1) = Easting;
    seat11_N(next_idx:next_idx+length(Northing)-1) = Northing;
    seat11_SMB_MC(next_idx:next_idx+length(SMB)-1) = SMB;
    seat11_yr(next_idx:next_idx+length(SMB_yr)-1) = SMB_yr;
    
end

% Combine SEAT2010 and SEAT2011 SMB data into combined arrays
seat_E = [seat10_E seat11_E];
seat_N = [seat10_N seat11_N];
seat_SMB_MC = [seat10_SMB_MC seat11_SMB_MC];
seat_yr = [seat10_yr seat11_yr];

% Find and remove empty SEAT indices (usually from excess preallocation)
keep_idx = find(~cellfun(@isempty, seat_SMB_MC));
seat_E = seat_E(keep_idx);
seat_N = seat_N(keep_idx);
seat_SMB_MC = seat_SMB_MC(keep_idx);
seat_yr = seat_yr(keep_idx);

% Calculate mean (and st. dev.) annual SMB from the MC simulations for 
% SEAT data
seat_SMB = cellfun(@(x) mean(x, 2), seat_SMB_MC, 'UniformOutput', 0);
seat_std = cellfun(@(x) std(x, [], 2), seat_SMB_MC, 'UniformOutput', 0);

% Load previously processed 2011 OIB snow radar accumulation results
OIB_files = dir(fullfile(data_path, 'IceBridge/SNO_radar/',...
    '2011/SMB_results/', wild));

oib_E = zeros(1, length(OIB_files)*2*(50*1000/25));
oib_N = oib_E;
oib_SMB_MC = cell(1, length(OIB_files)*2*(50*1000/25));
oib_yr = oib_SMB_MC;

for i = 1:length(OIB_files)
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Easting');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Northing');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB_yr');

    next_idx = sum(~cellfun(@isempty, oib_SMB_MC)) + 1;
    oib_E(next_idx:next_idx+length(Easting)-1) = Easting;
    oib_N(next_idx:next_idx+length(Northing)-1) = Northing;
    oib_SMB_MC(next_idx:next_idx+length(SMB)-1) = SMB;
    oib_yr(next_idx:next_idx+length(SMB_yr)-1) = SMB_yr;
end

keep_idx = find(~cellfun(@isempty, oib_SMB_MC));
oib_E = oib_E(keep_idx);
oib_N = oib_N(keep_idx);
oib_SMB_MC = oib_SMB_MC(keep_idx);
oib_yr = oib_yr(keep_idx);

oib_SMB = cellfun(@(x) mean(x, 2), oib_SMB_MC, 'UniformOutput', 0);
oib_std = cellfun(@(x) std(x, [], 2), oib_SMB_MC, 'UniformOutput', 0);

% Attempt to additionally load OIB elevation data, if available
try
    oib_elev = false(1, length(OIB_files)*2*(50*1000/25));
    for i=1:length(OIB_files)
        load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'elev');
        next_idx = sum(logical(oib_elev)) + 1;
        oib_elev(next_idx:next_idx+length(elev)-1) = elev;
    end
    keep_idx = find(logical(oib_elev));
    oib_elev = oib_elev(keep_idx);
catch
    disp('Flag: Missing elevation data')
end

%%

% Define starting and end year for regression, and create yr vector based
% on those values
yr_start = 1978;
yr_end = 2008;
year = (yr_end:-1:yr_start)';


%%

% Determine indices of SEAT data that covers the specified time period
seat_idx = cellfun(@(x) max(x)>=yr_end && min(x)<=yr_start, seat_yr);

% Extract SEAT data that covers specified time period
SEAT_E = seat_E(seat_idx);
SEAT_N = seat_N(seat_idx);
SEAT_SMB = seat_SMB(seat_idx);
SEAT_std = seat_std(seat_idx);
SEAT_yr = seat_yr(seat_idx);

% Determine the positions within the records of start and stop year for
% each SEAT trace
SEAT_start = cellfun(@(x) find(x==yr_start, 1), SEAT_yr, 'UniformOutput', 0);
SEAT_end = cellfun(@(x) find(x==yr_end, 1), SEAT_yr, 'UniformOutput', 0);

% Combine SEAT SMB data into matrix-convertible arrays (all cells are the
% same length and cover the same time period)
mSEAT_SMB = cellfun(@(x,y,z) x(y:z), SEAT_SMB, SEAT_end, SEAT_start, ...
    'UniformOutput', 0);
% mSEAT_SMB = cellfun(@(x,y,z) movmean(x(y:z),3), SEAT_SMB, SEAT_end, SEAT_start, ...
%     'UniformOutput', 0);

% Calculate the iteratively reweighted least squares regression for each
% SEAT trace time series and place in structure array 
[coeff, stats] = cellfun(@(x) robustfit(year, x), mSEAT_SMB, 'UniformOutput', 0);
SEAT_stats = struct();
SEAT_stats.b = cellfun(@(x) x(2), coeff);
SEAT_stats.se = cellfun(@(x) x.se(2), stats); 
SEAT_stats.p = cellfun(@(x) x.p(2), stats);

% [coeff, stats] = cellfun(@(x) polyfit(year, x, 1), mSEAT_SMB, 'UniformOutput', 0);
% SEAT_stats = struct();
% SEAT_stats.b = cellfun(@(x) x(1), coeff);


%%
oib_idx = cellfun(@(x) max(x)>=yr_end && min(x)<=yr_start, oib_yr);
OIB_E = oib_E(oib_idx);
OIB_N = oib_N(oib_idx);
OIB_SMB = oib_SMB(oib_idx);
OIB_std = oib_std(oib_idx);
OIB_yr = oib_yr(oib_idx);

OIB_start = cellfun(@(x) find(x==yr_start, 1), OIB_yr, 'UniformOutput', 0);
OIB_end = cellfun(@(x) find(x==yr_end, 1), OIB_yr, 'UniformOutput', 0);
mOIB_SMB = cellfun(@(x,y,z) x(y:z), OIB_SMB, OIB_end, OIB_start, ...
    'UniformOutput', 0);
% mOIB_SMB = cellfun(@(x,y,z) movmean(x(y:z),3), OIB_SMB, OIB_end, OIB_start, ...
%     'UniformOutput', 0);

[coeff, stats] = cellfun(@(x) robustfit(year, movmean(x,3)), mOIB_SMB, 'UniformOutput', 0);
OIB_stats = struct();
OIB_stats.b = cellfun(@(x) x(2), coeff);
OIB_stats.se = cellfun(@(x) x.se(2), stats); 
OIB_stats.p = cellfun(@(x) x.p(2), stats);

% [coeff, stats] = cellfun(@(x) polyfit(year, x, 1), mOIB_SMB, 'UniformOutput', 0);
% OIB_stats = struct();
% OIB_stats.b = cellfun(@(x) x(1), coeff);

%%

core_idx = 1:5;

cores_SMB = nan(length(year), length(cores.name(core_idx)));
cores_beta = zeros(1, length(cores.name(core_idx)));
cores_se = zeros(1, length(cores.name(core_idx)));
cores_pval = zeros(1, length(cores.name(core_idx)));
for k = 1:length(cores.name(core_idx))
    core_k = cores.(cores.name{k});
    core_start = find(core_k.SMB_yr<=yr_end, 1, 'first');
    core_end = find(core_k.SMB_yr>=yr_start, 1, 'last');
    SMB_mean = mean(core_k.SMB, 2);
    SMB_k = SMB_mean(core_start:core_end);
    cores_SMB(1:length(SMB_k),k) = SMB_k;
    
    [coeff, stats] = robustfit(year(1:length(SMB_k)), SMB_k);
    cores_beta(k) = coeff(2);
    cores_se(k) = stats.se(2);
    cores_pval(k) = stats.p(2);

% [coeff, stats] = polyfit(year(1:length(SMB_k)), SMB_k, 1);
%     cores_beta(k) = coeff(1);

end

%% Main map of mean SMB

labels = strrep(cores.name(core_idx), '_', '-');
basins = shaperead(strcat(data_path, ...
    'DEMs/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'));
Easting_lims = [min([min(cores.Easting) min(SEAT_E) min(OIB_E)]) - 5000 ...
    max([max(cores.Easting) max(SEAT_E) max(OIB_E)]) + 5000];
Northing_lims = [min([min(cores.Northing) min(SEAT_N) min(OIB_N)]) - 5000 ...
    max([max(cores.Northing) max(SEAT_N) max(OIB_N)]) + 5000];

[Arth_E, Arth_N, Arth_accum] = accumulation_data(Easting_lims, Northing_lims, 'xy');


map_SMB = figure('Position', [10 10 1400 800]);
h0 = image(Arth_E(1,:), (Arth_N(:,1))', Arth_accum, 'CDataMapping', 'scaled');
set(gca, 'Ydir', 'normal')
hold on
h1 = mapshow(basins, 'FaceAlpha', 0);
h2 = scatter(SEAT_E, SEAT_N, 50, mean(cell2mat(mSEAT_SMB)), 'filled');
h3 = scatter(OIB_E, OIB_N, 50, mean(cell2mat(mOIB_SMB)), 'filled');
h4 = scatter(cores.Easting(core_idx), cores.Northing(core_idx), 125, ...
    nanmean(cores_SMB)', 'filled', 'MarkerEdgeColor', 'k');
text(cores.Easting(core_idx), cores.Northing(core_idx), ...
    strcat(labels, '\rightarrow'), 'FontSize', 13, ...
    'Interpreter', 'tex', 'HorizontalAlignment', 'right');
c0 = colorbar;
caxis([150 450])
c0.Label.String = ['Mean annual SMB ' num2str(yr_start) '-' ...
    num2str(yr_end) ' (mm/a)'];
c0.Label.FontSize = 18;
graticuleps(-81:0.5:-77,-125:2:-105, 'c')
xlim(Easting_lims)
ylim(Northing_lims)
scalebarps
box on
mapzoomps('ne', 'insetsize', 0.30)
% legend([h0 h3 h4], 'Arthern mean SMB', 'SEAT core mean SMB', ...
%     'SEAT radar mean SMB', 'Location', 'northwest')
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 18)
% title('SEAT mean annual SMB')
hold off

map1_name = 'SMB_mean_map';
export_fig(map_SMB, fullfile(output_dir, map1_name), '-png');
close(map_SMB)


%% Main map of SMB trends

% Add addon to generate custom color scale to path
addon_folder = fullfile(addon_path, 'b2r');
addpath(genpath(addon_folder))

map_trend = figure('Position', [10 10 1400 800]);
% title('SEAT radar SMB trends')
hold on
h1 = mapshow(basins, 'FaceAlpha', 0);
h2 = scatter(SEAT_E, SEAT_N, 50, SEAT_stats.b, 'filled');
% h2 = scatter(SEAT_E(SEAT_stats.p<=0.05), SEAT_N(SEAT_stats.p<=0.05), 50, ...
%     SEAT_stats.b(SEAT_stats.p<=0.05), 'filled');
h3 = scatter(SEAT_E(SEAT_stats.p<=0.05), SEAT_N(SEAT_stats.p<=0.05), 3, ...
    'y', 'filled', 'MarkerFaceAlpha', 0.25, ...
    'MarkerEdgeAlpha', 0.25);
h4 = scatter(OIB_E, OIB_N, 50, OIB_stats.b, 'filled');
% h4 = scatter(OIB_E(OIB_stats.p<=0.05), OIB_N(OIB_stats.p<=0.05), 50, ...
%     OIB_stats.b(OIB_stats.p<=0.05), 'filled');
h5 = scatter(OIB_E(OIB_stats.p<=0.05), OIB_N(OIB_stats.p<=0.05), 3, ...
    'y', 'filled', 'MarkerFaceAlpha', 0.25, ...
    'MarkerEdgeAlpha', 0.25);
h6 = scatter(cores.Easting(core_idx), cores.Northing(core_idx), 125, ...
    cores_beta, 'filled', 'MarkerEdgeColor', 'k');
text(cores.Easting, cores.Northing, strcat(labels, '\rightarrow'), ...
    'FontSize', 13, 'Interpreter', 'tex', 'HorizontalAlignment', 'right');
colormap(b2r(-4, 2))
c0 = colorbar;
c0.Label.String = ['Linear trend in annual SMB ' num2str(yr_start) '-' ...
    num2str(yr_end) ' (mm/a)'];
c0.Label.FontSize = 18;
graticuleps(-81:0.5:-77,-125:2:-105, 'c')
xlim(Easting_lims)
ylim(Northing_lims)
scalebarps
box on
mapzoomps('ne', 'insetsize', 0.30)
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 18)
hold off

%%

