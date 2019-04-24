% This script generates the results figures and values for comparisons
% between OIB and SEAT radar in the TGARS paper (Keeler et al, 2019)

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
                
            case 'laptop'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = fullfile('C:/Users/durba/', ...
                    'OneDrive - University of Utah/Documents/MATLAB/Addons/');
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
SEAT_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
    'SEAT2010Kuband/SEAT10_4toSEAT10_6/SMB_results/', wild));

% Preallocate arrays of sufficient size for data
SEAT_E = zeros(1, length(SEAT_files)*2*(50*1000/25));
SEAT_N = SEAT_E;
SEAT_SMB_MC = cell(1, length(SEAT_files)*2*(50*1000/25));
SEAT_yr = SEAT_SMB_MC;

for i = 1:length(SEAT_files)
    
    % Load relevent data from current data file
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'Easting');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'Northing');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'SMB');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'SMB_yr');
    
    % Find position of last data entered into preallocated arrays
    next_idx = sum(~cellfun(@isempty, SEAT_SMB_MC)) + 1;
    
    % Fill current iteration data into preallocated arrays
    SEAT_E(next_idx:next_idx+length(Easting)-1) = Easting;
    SEAT_N(next_idx:next_idx+length(Northing)-1) = Northing;
    SEAT_SMB_MC(next_idx:next_idx+length(SMB)-1) = SMB;
    SEAT_yr(next_idx:next_idx+length(SMB_yr)-1) = SMB_yr;
    
end

% Find and remove empty SEAT indices (usually from excess preallocation)
keep_idx = find(~cellfun(@isempty, SEAT_SMB_MC));
SEAT_E = SEAT_E(keep_idx);
SEAT_N = SEAT_N(keep_idx);
SEAT_SMB_MC = SEAT_SMB_MC(keep_idx);
SEAT_yr = SEAT_yr(keep_idx);

% Calculate mean (and st. dev.) annual SMB from the MC simulations for
% SEAT data
seat_SMB = cellfun(@(x) mean(x, 2), SEAT_SMB_MC, 'UniformOutput', 0);
seat_std = cellfun(@(x) std(x, [], 2), SEAT_SMB_MC, 'UniformOutput', 0);

% Load previously processed 2011 OIB snow radar accumulation results
OIB_files = dir(fullfile(data_path, 'IceBridge/SEAT10_4to10_6/',...
    '2011_SNO/SMB_results/', wild));

OIB_E = zeros(1, length(OIB_files)*2*(50*1000/25));
OIB_N = OIB_E;
OIB_SMB_MC = cell(1, length(OIB_files)*2*(50*1000/25));
OIB_yr = OIB_SMB_MC;

for i = 1:length(OIB_files)
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Easting');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Northing');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB_yr');
    
    next_idx = sum(~cellfun(@isempty, OIB_SMB_MC)) + 1;
    OIB_E(next_idx:next_idx+length(Easting)-1) = Easting;
    OIB_N(next_idx:next_idx+length(Northing)-1) = Northing;
    OIB_SMB_MC(next_idx:next_idx+length(SMB)-1) = SMB;
    OIB_yr(next_idx:next_idx+length(SMB_yr)-1) = SMB_yr;
end

keep_idx = find(~cellfun(@isempty, OIB_SMB_MC));
OIB_E = OIB_E(keep_idx);
OIB_N = OIB_N(keep_idx);
OIB_SMB_MC = OIB_SMB_MC(keep_idx);
OIB_yr = OIB_yr(keep_idx);

oib_SMB = cellfun(@(x) mean(x, 2), OIB_SMB_MC, 'UniformOutput', 0);
oib_std = cellfun(@(x) std(x, [], 2), OIB_SMB_MC, 'UniformOutput', 0);

% Attempt to additionally load OIB elevation data, if available
try
    OIB_elev = false(1, length(OIB_files)*2*(50*1000/25));
    for i=1:length(OIB_files)
        load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'elev');
        next_idx = sum(logical(OIB_elev)) + 1;
        OIB_elev(next_idx:next_idx+length(elev)-1) = elev;
    end
    keep_idx = find(logical(OIB_elev));
    OIB_elev = OIB_elev(keep_idx);
catch
    disp('Flag: Missing elevation data')
end

%% Study site map

% Get subset of cores used in study
core_idx = 3:5;

% Generate labels
labels = strrep(cores.name(core_idx), '_', '-');

% Get basin boundaries, and generate map limits based on included data
% basins = shaperead(strcat(data_path, ...
%     'DEMs/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'));
Easting_lims = [min([min(cores.Easting(core_idx)) min(SEAT_E) min(OIB_E)]) ...
    - 15000 max([max(cores.Easting(core_idx)) max(SEAT_E) max(OIB_E)]) + 15000];
Northing_lims = [min([min(cores.Northing(core_idx)) min(SEAT_N) min(OIB_N)])...
    - 30000 max([max(cores.Northing(core_idx)) max(SEAT_N) max(OIB_N)]) + 30000];

% Generate bulk accumulation rate estimates from Arthern et al (2006)
[Arth_E, Arth_N, Arth_accum] = accumulation_data(Easting_lims, Northing_lims, 'xy');

% Get elevation data
[elev_E, elev_N, elev] = cryosat2_data(Easting_lims, Northing_lims, 'xy');

% Plot of the study site
fig = figure('Position', [0 0 1400 800]);
h0 = image(Arth_E(1,:), (Arth_N(:,1))', Arth_accum, 'CDataMapping', 'scaled');
set(gca, 'Ydir', 'normal')
hold on
h1 = contour(elev_E(1,:), elev_N(:,1)', elev, 'k', 'showtext', 'on');
set(gca,'clim',[min(Arth_accum(:)) max(Arth_accum(:))]);
% h2 = plot(SEAT_E, SEAT_N, 'r.');
h3 = plot(OIB_E, OIB_N, 'm.');
h4 = scatter(cores.Easting(core_idx), cores.Northing(core_idx), 50, 'b', ...
    'filled', 'MarkerEdgeColor', 'k');
text(cores.Easting(core_idx), cores.Northing(core_idx), labels, 'white', ...
    'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
c0 = colorbar;
c0.Label.String = 'Mean annual accumulation (mm/a)';
c0.Label.FontSize = 14;
xlim(Easting_lims)
ylim(Northing_lims)
scalebarps
box on
mapzoomps('ne', 'insetsize', 0.30)
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 12)
hold off
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 4], ...
    'PaperUnits', 'Inches', 'PaperSize', [8, 4])

% % Export figure
% fig_nm = 'site-map';
% export_fig(fig, fullfile(output_dir, fig_nm), '-pdf', '-q101', '-cmyk', '-a1')
% close(fig)


%% SEAT and OIB SMB bias tests

SEAT_SMB = cellfun(@(x) mean(x, 2), SEAT_SMB_MC, 'UniformOutput', 0);
OIB_SMB = cellfun(@(x) mean(x, 2), OIB_SMB_MC, 'UniformOutput', 0);

near_tmp = knnsearch([SEAT_E' SEAT_N'], [OIB_E' OIB_N']);
OIB_near = 1:length(OIB_E);
SEAT_near = near_tmp;
D_near = diag(pdist2([OIB_E(OIB_near)' OIB_N(OIB_near)'], ...
    [SEAT_E(SEAT_near)' SEAT_N(SEAT_near)']));
D_idx = D_near<=25;
OIB_near = OIB_near(D_idx);
SEAT_near = SEAT_near(D_idx);

yr_start = cellfun(@(x,y) min([x(1) y(1)]), ...
    OIB_yr(OIB_near), SEAT_yr(SEAT_near));
yr_end = cellfun(@(x,y) max([x(end), y(end)]), ...
    OIB_yr(OIB_near), SEAT_yr(SEAT_near));

bias_yr = cell(1, length(yr_start));
for j = 1:length(yr_start)
    bias_yr{j} = (yr_start(j):-1:yr_end(j))';
end

SEAT_start = cellfun(@(x,y) find(x==y(1), 1), SEAT_yr(SEAT_near), bias_yr);
SEAT_end = cellfun(@(x,y) find(x==y(end), 1), SEAT_yr(SEAT_near), bias_yr);
OIB_start = cellfun(@(x,y) find(x==y(1), 1), OIB_yr(OIB_near), bias_yr);
OIB_end = cellfun(@(x,y) find(x==y(end), 1), OIB_yr(OIB_near), bias_yr);

SEAT_SMB_mu = cell(1,length(bias_yr));
SMB_res = cell(1,length(bias_yr));
for j=1:length(bias_yr)
    SEAT_SMB_j = SEAT_SMB{SEAT_near(j)}(SEAT_start(j):SEAT_end(j));
%     SEATbias_mean{j} = mean(SEATbias_SMB);
    SEAT_SMB_mu{j} = SEAT_SMB_j;
    OIB_SMB_j = OIB_SMB{OIB_near(j)}(OIB_start(j):OIB_end(j));
    SMB_res{j} = (SEAT_SMB_j - OIB_SMB_j);
end


SEAT_SMB_mu = cell(1,length(bias_yr));
SMB_res = cell(1,length(bias_yr));
for j=1:length(bias_yr)
    SEAT_SMB_j = SEAT_SMB{SEAT_near(j)}(SEAT_start(j):SEAT_end(j));
    SEAT_SMB_mu{j} = SEAT_SMB_j;
    OIB_SMB_j = OIB_SMB{OIB_near(j)}(OIB_start(j):OIB_end(j));
    SMB_res{j} = (SEAT_SMB_j - OIB_SMB_j);
end


% SMB bias stats for each individual year for each mean trace within the
% threshold distance (% difference of mean SMB)
bias_dist = cellfun(@(x,y) x./y, SMB_res, SEAT_SMB_mu, 'UniformOutput', 0);
bias_dist = vertcat(bias_dist{:});
MBE = median(bias_dist);
MBE_mu = mean(bias_dist);
RMSE = sqrt(mean(bias_dist.^2));
bias_stats.radar = table(MBE, MBE_mu, RMSE, 'VariableNames', ...
    {'MedianMBE', 'MeanMBE', 'RMSE'}, 'RowNames', {'SEAT-OIB (% bias)'});

% bias_std = std(bias_dist);
% bias_SEM = bias_std/sqrt(length(bias_dist));
% Ts = tinv(0.975, length(bias_dist)-1);
% bias_MoE = Ts*bias_SEM;
% figure
% histogram(bias_dist, 100)
% bias_stats.radar = table(bias_med, bias_mu, bias_MoE, bias_std, 'VariableNames', ...
%     {'Median', 'Mean','MarginOfError','StdDev'}, 'RowNames', {'SEAT-OIB (% bias)'});





SEAT_var = cellfun(@(x) var(x,[],2)./(mean(x,2)).^2, SEAT_SMB_MC(SEAT_near), ...
    'UniformOutput', false);
SEAT_std = sqrt(median(vertcat(SEAT_var{:})));
% std(vertcat(SEAT_std{:}));
OIB_var = cellfun(@(x) var(x,[],2)./(mean(x,2)).^2, OIB_SMB_MC(OIB_near), ...
    'UniformOutput', false);
OIB_std = sqrt(median(vertcat(OIB_var{:})));

clip30 = cellfun(@(x,y) length(x)>=30 && length(y)>=30, SEAT_var, OIB_var); 
SEAT_std30 = sqrt(mean(cell2mat(cellfun(@(x) x(1:30), SEAT_var(clip30), ...
    'UniformOutput', false)), 2));
OIB_std30 = sqrt(mean(cell2mat(cellfun(@(x) x(1:30), OIB_var(clip30), ...
    'UniformOutput', false)), 2));
