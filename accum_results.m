% File to produce final estimates and results from the SEAT2010-4 to
% SEAT2010-6 radar line

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

%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);
Ndraw = 100;

% Load previously processed SEAT2010 accumulation data
wild = '*.mat';
SEAT10_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
    'SEAT2010Kuband/RawSEAT2010/SMB_results/', wild));

seat10_E = zeros(1, length(SEAT10_files)*2*(50*1000/25));
seat10_N = seat10_E;
seat10_SMB_MC = cell(1, length(SEAT10_files)*2*(50*1000/25));
seat10_yr = seat10_SMB_MC;
for i = 1:length(SEAT10_files)
    load(fullfile(SEAT10_files(i).folder, SEAT10_files(i).name), 'Easting');
    load(fullfile(SEAT10_files(i).folder, SEAT10_files(i).name), 'Northing');
    load(fullfile(SEAT10_files(i).folder, SEAT10_files(i).name), 'SMB');
    load(fullfile(SEAT10_files(i).folder, SEAT10_files(i).name), 'SMB_yr');
    
    next_idx = sum(~cellfun(@isempty, seat10_SMB_MC)) + 1;
    seat10_E(next_idx:next_idx+length(Easting)-1) = Easting;
    seat10_N(next_idx:next_idx+length(Northing)-1) = Northing;
    seat10_SMB_MC(next_idx:next_idx+length(SMB)-1) = SMB;
    seat10_yr(next_idx:next_idx+length(SMB_yr)-1) = SMB_yr;
    
end

SEAT11_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
    'SEAT2011Kuband/RawSEAT2011/SMB_results/', wild));

seat11_E = zeros(1, length(SEAT11_files)*2*(50*1000/25));
seat11_N = seat11_E;
seat11_SMB_MC = cell(1, length(SEAT11_files)*2*(50*1000/25));
seat11_yr = seat10_SMB_MC;
for i = 1:length(SEAT11_files)
    load(fullfile(SEAT11_files(i).folder, SEAT11_files(i).name), 'Easting');
    load(fullfile(SEAT11_files(i).folder, SEAT11_files(i).name), 'Northing');
    load(fullfile(SEAT11_files(i).folder, SEAT11_files(i).name), 'SMB');
    load(fullfile(SEAT11_files(i).folder, SEAT11_files(i).name), 'SMB_yr');
    
    next_idx = sum(~cellfun(@isempty, seat11_SMB_MC)) + 1;
    seat11_E(next_idx:next_idx+length(Easting)-1) = Easting;
    seat11_N(next_idx:next_idx+length(Northing)-1) = Northing;
    seat11_SMB_MC(next_idx:next_idx+length(SMB)-1) = SMB;
    seat11_yr(next_idx:next_idx+length(SMB_yr)-1) = SMB_yr;
    
end

seat_E = [seat10_E seat11_E];
seat_N = [seat10_N seat11_N];
seat_SMB_MC = [seat10_SMB_MC seat11_SMB_MC];
seat_yr = [seat10_yr seat11_yr];


keep_idx = find(~cellfun(@isempty, seat_SMB_MC));
seat_E = seat_E(keep_idx);
seat_N = seat_N(keep_idx);
seat_SMB_MC = seat_SMB_MC(keep_idx);
seat_yr = seat_yr(keep_idx);

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

yr_start = 1979;
yr_end = 2009;
year = (yr_end:-1:yr_start)';


%%
seat_idx = cellfun(@(x) max(x)>=yr_end && min(x)<=yr_start, seat_yr);
SEAT_E = seat_E(seat_idx);
SEAT_N = seat_N(seat_idx);
SEAT_SMB = seat_SMB(seat_idx);
SEAT_std = seat_std(seat_idx);
SEAT_yr = seat_yr(seat_idx);

SEAT_start = cellfun(@(x) find(x==yr_start, 1), SEAT_yr, 'UniformOutput', 0);
SEAT_end = cellfun(@(x) find(x==yr_end, 1), SEAT_yr, 'UniformOutput', 0);
% mSEAT_SMB = cellfun(@(x,y,z) x(y:z), SEAT_SMB, SEAT_end, SEAT_start, ...
%     'UniformOutput', 0);
mSEAT_SMB = cellfun(@(x,y,z) movmean(x(y:z),3), SEAT_SMB, SEAT_end, SEAT_start, ...
    'UniformOutput', 0);

% [coeff, stats] = cellfun(@(x) robustfit(year, x), mSEAT_SMB, 'UniformOutput', 0);
% SEAT_stats = struct();
% SEAT_stats.b = cellfun(@(x) x(2), coeff);
% SEAT_stats.se = cellfun(@(x) x.se(2), stats); 
% SEAT_stats.p = cellfun(@(x) x.p(2), stats);

[coeff, stats] = cellfun(@(x) polyfit(year, x, 1), mSEAT_SMB, 'UniformOutput', 0);
SEAT_stats = struct();
SEAT_stats.b = cellfun(@(x) x(1), coeff);
% SEAT_stats.se = cellfun(@(x) x.se(2), stats); 
% SEAT_stats.p = cellfun(@(x) x.p(2), stats);

%%
oib_idx = cellfun(@(x) max(x)>=yr_end && min(x)<=yr_start, oib_yr);
OIB_E = oib_E(oib_idx);
OIB_N = oib_N(oib_idx);
OIB_SMB = oib_SMB(oib_idx);
OIB_std = oib_std(oib_idx);
OIB_yr = oib_yr(oib_idx);

OIB_start = cellfun(@(x) find(x==yr_start, 1), OIB_yr, 'UniformOutput', 0);
OIB_end = cellfun(@(x) find(x==yr_end, 1), OIB_yr, 'UniformOutput', 0);
% mOIB_SMB = cellfun(@(x,y,z) x(y:z), OIB_SMB, OIB_end, OIB_start, ...
%     'UniformOutput', 0);
mOIB_SMB = cellfun(@(x,y,z) movmean(x(y:z),3), OIB_SMB, OIB_end, OIB_start, ...
    'UniformOutput', 0);

% [coeff, stats] = cellfun(@(x) robustfit(year, movmean(x,3)), mOIB_SMB, 'UniformOutput', 0);
% OIB_stats = struct();
% OIB_stats.b = cellfun(@(x) x(2), coeff);
% OIB_stats.se = cellfun(@(x) x.se(2), stats); 
% OIB_stats.p = cellfun(@(x) x.p(2), stats);

[coeff, stats] = cellfun(@(x) polyfit(year, x, 1), mOIB_SMB, 'UniformOutput', 0);
OIB_stats = struct();
OIB_stats.b = cellfun(@(x) x(1), coeff);
% OIB_stats.se = cellfun(@(x) x.se(2), stats); 
% OIB_stats.p = cellfun(@(x) x.p(2), stats);

%%

cores_SMB = nan(length(year), length(cores.name));
cores_beta = zeros(1, length(cores.name));
cores_se = zeros(1, length(cores.name));
cores_pval = zeros(1, length(cores.name));
for k = 1:length(cores.name)
    core_k = cores.(cores.name{k});
    core_start = find(core_k.SMB_yr<=yr_end, 1, 'first');
    core_end = find(core_k.SMB_yr>=yr_start, 1, 'last');
    SMB_mean = mean(core_k.SMB, 2);
    SMB_k = SMB_mean(core_start:core_end);
    cores_SMB(1:length(SMB_k),k) = SMB_k;
    
%     [coeff, stats] = robustfit(year(1:length(SMB_k)), SMB_k);
%     cores_beta(k) = coeff(2);
%     cores_se(k) = stats.se(2);
%     cores_pval(k) = stats.p(2);

[coeff, stats] = polyfit(year(1:length(SMB_k)), SMB_k, 1);
    cores_beta(k) = coeff(1);
%     cores_se(k) = stats.se(2);
%     cores_pval(k) = stats.p(2);
end


%% Index map and study site

labels = strrep(cores.name, '_', '-');
basins = shaperead(strcat(data_path, ...
    'DEMs/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'));
Easting_lims = [min([min(cores.Easting) min(SEAT_E) min(OIB_E)]) - 5000 ...
    max([max(cores.Easting) max(SEAT_E) max(OIB_E)]) + 5000];
Northing_lims = [min([min(cores.Northing) min(SEAT_N) min(OIB_N)]) - 5000 ...
    max([max(cores.Northing) max(SEAT_N) max(OIB_N)]) + 5000];

[Arth_E, Arth_N, Arth_accum] = accumulation_data(Easting_lims, Northing_lims, 'xy');


figure('Position', [10 10 1400 800])
h0 = image(Arth_E(1,:), (Arth_N(:,1))', Arth_accum, 'CDataMapping', 'scaled');
set(gca, 'Ydir', 'normal')
hold on
h1 = mapshow(basins, 'FaceAlpha', 0);
h2 = scatter(SEAT_E, SEAT_N, 25, mean(cell2mat(mSEAT_SMB)), 'filled');
% h3 = scatter(OIB_E, OIB_N, 25, mean(cell2mat(mOIB_SMB)), 'filled');
h4 = scatter(cores.Easting, cores.Northing, 100, nanmean(cores_SMB)', 'filled');
text(cores.Easting, cores.Northing, strcat('\leftarrow', labels), ...
    'FontSize', 13, 'Interpreter', 'tex');
c0 = colorbar;
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
title('SEAT mean annual SMB')
hold off


figure('Position', [10 10 1400 800])
title('SEAT radar SMB trends')
hold on
h1 = mapshow(basins, 'FaceAlpha', 0);
h2 = scatter(SEAT_E, SEAT_N, 25, SEAT_stats.b, 'filled');
% h3 = scatter(OIB_E, OIB_N, 25, OIB_stats.b, 'filled');
h4 = scatter(cores.Easting, cores.Northing, 100, cores_beta, 'filled');
text(cores.Easting, cores.Northing, strcat('\leftarrow', labels), ...
    'FontSize', 15, 'Interpreter', 'tex');
c0 = colorbar;
c0.Label.String = ['Linear trend in annual SMB ' num2str(yr_start) '-' ...
    num2str(yr_end) ' (mm/a)'];
c0.Label.FontSize = 18;
graticuleps(-81:0.5:-77,-125:2:-105, 'c')
xlim(Easting_lims)
ylim(Northing_lims)
caxis([-7 0])
scalebarps
box on
mapzoomps('ne', 'insetsize', 0.30)
% legend([h0 h3 h4], 'Arthern mean SMB', 'SEAT core mean SMB', ...
%     'SEAT radar mean SMB', 'Location', 'northwest')
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 18)
hold off

%%

figure('Position', [10 10 1400 800])
h0 = image(Arth_E(1,:), (Arth_N(:,1))', Arth_accum, 'CDataMapping', 'scaled');
set(gca, 'Ydir', 'normal')
hold on
h1 = mapshow(basins, 'FaceAlpha', 0);
% h2 = scatter(SEAT_E, SEAT_N, 25, mean(cell2mat(mSEAT_SMB)), 'filled');
h3 = scatter(OIB_E, OIB_N, 25, mean(cell2mat(mOIB_SMB)), 'filled');
h4 = scatter(cores.Easting, cores.Northing, 100, nanmean(cores_SMB)', 'filled');
text(cores.Easting, cores.Northing, strcat('\leftarrow', labels), ...
    'FontSize', 13, 'Interpreter', 'tex');
c0 = colorbar;
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
title('OIB mean annual SMB')
hold off


figure('Position', [10 10 1400 800])
title('OIB radar SMB trends')
hold on
h1 = mapshow(basins, 'FaceAlpha', 0);
% h2 = scatter(SEAT_E, SEAT_N, 25, SEAT_beta, 'filled');
h3 = scatter(OIB_E, OIB_N, 25, OIB_stats.b, 'filled');
h4 = scatter(cores.Easting, cores.Northing, 100, cores_beta, 'filled');
text(cores.Easting, cores.Northing, strcat('\leftarrow', labels), ...
    'FontSize', 15, 'Interpreter', 'tex');
c0 = colorbar;
c0.Label.String = ['Linear trend in annual SMB ' num2str(yr_start) '-' ...
    num2str(yr_end) ' (mm/a)'];
c0.Label.FontSize = 18;
graticuleps(-81:0.5:-77,-125:2:-105, 'c')
xlim(Easting_lims)
ylim(Northing_lims)
caxis([-7 0])
scalebarps
box on
mapzoomps('ne', 'insetsize', 0.30)
% legend([h0 h3 h4], 'Arthern mean SMB', 'SEAT core mean SMB', ...
%     'SEAT radar mean SMB', 'Location', 'northwest')
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 18)
hold off






%%

SEAT_set = SEAT_E >= -1.05E6 & SEAT_E <= -1.03E6 & SEAT_N <= -4.638E5;
OIB_set = OIB_E >= -1.05E6 & OIB_E <= -1.03E6 & OIB_N <= -4.638E5;
SMB_S = SEAT_SMB(SEAT_set);
SMB_O = fliplr(OIB_SMB(OIB_set));
std_S = SEAT_std(SEAT_set);
std_O = fliplr(OIB_std(OIB_set));
yr_S = SEAT_yr(SEAT_set);
yr_O = fliplr(OIB_yr(OIB_set));

% figure
% hold on
% scatter(SEAT_E(SEAT_set), SEAT_N(SEAT_set), 25, ...
%     SEAT_beta(SEAT_set), 'filled')
% scatter(OIB_E(OIB_set), OIB_N(OIB_set), 25, OIB_beta(OIB_set), 'filled')
% hold off

i = randi(length(SMB_O));
figure
hold on
plot(yr_S{i}, movmean(SMB_S{i},5), 'r', 'LineWidth', 2)
plot(yr_S{i}, movmean(SMB_S{i},5) + movmean(std_S{i},5), 'r--')
plot(yr_S{i}, movmean(SMB_S{i},5) - movmean(std_S{i},5), 'r--')
plot(yr_O{i}, movmean(SMB_O{i},5), 'm', 'LineWidth', 2)
plot(yr_O{i}, movmean(SMB_O{i},5) + movmean(std_O{i},5), 'm--')
plot(yr_O{i}, movmean(SMB_O{i},5) - movmean(std_O{i},5), 'm--')
hold off

%%





%%

% %% Core site comparisons
% 
% inputs = {'SEAT10_4', 'SEAT10_5', 'SEAT10_6'};
% out_dir = fullfile('/media/durbank/CUmind/proposal/Figures/');
% trends = struct(inputs{1}, [], inputs{2}, [], inputs{3}, []);
% p_vals = struct(inputs{1}, [], inputs{2}, [], inputs{3}, []);
% 
% for i = 1:length(inputs)
%     
%     % Name of ith core site to perform analysis/stats on
%     name = inputs{i};
%     file = fullfile(data_path, 'radar/SEAT_Traverses/results_data', ...
%         strcat('grid', name, '.mat'));
%     radar_i = load(file);
%     
%     % Find OIB data within 5 km of the ith core and the nearest trace to
%     % the core
%     dist_OIB = pdist2([cores.(name).Easting  cores.(name).Northing], ...
%         [radar_OIB.Easting', radar_OIB.Northing']);
%     trace_idx = 1:length(radar_OIB.SMB);
%     OIB_idx = trace_idx(dist_OIB<=5000);
%     [~, OIB_near] = min(dist_OIB);
%     
% %     % Find full SEAT data within 7.5 km of the ith core, and the nearest
% %     % trace to the core
% %     dist_SEAT = pdist2([cores.(name).Easting  cores.(name).Northing], ...
% %         [radar_SEAT.Easting', radar_SEAT.Northing']);
% %     trace_idx = 1:length(radar_SEAT.SMB);
% %     SEAT_idx = trace_idx(dist_SEAT<=7500);
% %     [~, SEAT_near] = min(dist_SEAT);
%     
%     % Find the SEAT core site radar trace nearest the ith core
%     dist_SEATi = pdist2([cores.(name).Easting  cores.(name).Northing], ...
%         [radar_i.Easting', radar_i.Northing']);
%     [~, SEATi_near] = min(dist_SEATi);
%     
%     % Find manual layer count data within 5 km of the ith core, along with
%     % the trace nearest to the core
%     dist_man = pdist2([cores.(name).Easting  cores.(name).Northing], ...
%         [radar_man.Easting', radar_man.Northing']);
%     trace_idx = 1:length(radar_man.Easting);
%     man_idx = trace_idx(dist_man<=5000);
%     if isempty(man_idx)
%         man_near = man_idx;
%     else
%         [~, man_near] = min(dist_man);
%     end
%     man_cutoff = median(sum(~logical(isnan(radar_man.ages(:,man_idx)))));
%     
%     % Define naming conventions for this iteration
%     text_name = strrep(name, '_', '-');
%     
% %     % Plot of the OIB radargram within range of core, along with the
% %     % position of the ith core
% %     f1 = figure('Position', [200 200 1000 550]);
% %     imagesc(radar_OIB.Easting(OIB_idx), radar_OIB.depth, radar_OIB.data_smooth(:,OIB_idx), [-2 2])
% %     colorbar
% %     xlabel('Distance along profile (m)')
% %     ylabel('Depth (m)')
% %     hold on
% %     plot([cores.(name).Easting cores.(name).Easting], ...
% %         [radar_OIB.depth(1) radar_OIB.depth(end)], 'b', 'LineWidth', 2)
% %     set(gca, 'Ydir', 'reverse', 'FontSize', 14)
% %     title(strcat(text_name, ' radargram'));
% %     hold off
% %     f1_name = strcat(name, '_radargram');
% %     %     export_fig(f1, strcat(out_dir, f1_name), '-png');
% %     %     close(f1)
%     
% %     % Plot of age-depth estimates from OIB, SEAT core site, manual counts,
% %     % and core data (all realizations)
% %     f2 = figure('Position', [200 200 750 550]);
% %     hold on
% %     for n = OIB_idx
% %         h0 = plot(radar_OIB.depth, median(radar_OIB.ages(:,n,:), 3), 'm', 'LineWidth', 0.5);
% %         h0.Color(4) = 0.02;
% %     end
% %     for n = 1:length(radar_i.SMB)
% %         h0 = plot(radar_i.depth, median(radar_i.ages(:,n,:), 3), 'r', 'LineWidth', 0.5);
% %         h0.Color(4) = 0.02;
% %     end
% % %     for n = SEAT_idx
% % %         h0 = plot(radar_SEAT.depth, median(radar_SEAT.ages(:,n,:), 3), 'r', 'LineWidth', 0.5);
% % %         h0.Color(4) = 0.02;
% % %     end
% %     for n = man_idx
% %         h0 = plot(radar_man.depth, radar_man.ages(:,n), 'k', 'LineWidth', 0.5);
% %         h0.Color(4) = 0.02;
% %     end
% %     h1 = plot(radar_i.depth, median(median(radar_i.ages, 3), 2),...
% %         'r--', 'LineWidth', 1);
% %     h2 = plot(radar_OIB.depth, median(median(radar_OIB.ages(:,OIB_idx,:), 3), 2),...
% %         'm--', 'LineWidth', 1);
% %     h3 = plot(cores.(name).depth, mean(cores.(name).ages, 2), 'b', 'LineWidth', 2);
% %     plot(cores.(name).depth, mean(cores.(name).ages, 2) + ...
% %         2*std(cores.(name).ages, [], 2), 'b--')
% %     plot(cores.(name).depth, mean(cores.(name).ages, 2) - ...
% %         2*std(cores.(name).ages, [], 2), 'b--')
% %     try
% %         h4 = plot(radar_man.depth(1:man_cutoff), ...
% %             median(radar_man.ages(1:man_cutoff,man_idx), 2, 'omitnan'), 'k');
% %         legend([h1 h2 h3 h4],'SEAT traces', 'OIB traces', 'Core', 'SEAT manual')
% %     catch
% %         legend([h1 h2 h3],'SEAT traces', 'OIB traces', 'Core')
% %     end
% %     title(strcat(text_name, ' age-depth scale'))
% %     ylabel('Calendar years')
% %     xlabel('Depth (m)')
% %     hold off
% %     f2_name = strcat(name, '_age');
% %     %     export_fig(f1, strcat(out_dir, f2_name), '-png');
% %     %     close(f2)
%     
%     % Plot of age-depth estimates from OIB, SEAT core site, manual counts,
%     % and core data (only nearest trace)
%     f1 = figure('Position', [200 200 700 700]);
%     hold on
%         for n = 1:Ndraw
%             h0 = plot(radar_OIB.depth, radar_OIB.ages(:,OIB_near,n), 'm', 'LineWidth', 0.5);
%             h0.Color(4) = 0.05;
%         end
%         for n = 1:Ndraw
%             h0 = plot(radar_i.depth, radar_i.ages(:,SEATi_near,n), 'r', 'LineWidth', 0.5);
%             h0.Color(4) = 0.05;
%         end
%         for n = 1:Ndraw
%             h0 = plot(cores.(name).depth, cores.(name).ages(:,n), 'b', 'LineWidth', 0.5);
%             h0.Color(4) = 0.05;
%         end
%     h1 = plot(radar_i.depth, median(radar_i.ages(:,SEATi_near,:), 3),...
%         'r', 'LineWidth', 2);
% %     plot(radar_i.depth, median(radar_i.ages(:,SEATi_near,:), 3)...
% %         + 2*std(squeeze(radar_i.ages(:,SEATi_near,:)), [], 2), 'r--')
% %     plot(radar_i.depth, median(radar_i.ages(:,SEATi_near,:), 3)...
% %         - 2*std(squeeze(radar_i.ages(:,SEATi_near,:)), [], 2), 'r--')
%     h2 = plot(radar_OIB.depth, median(radar_OIB.ages(:,OIB_near,:), 3),...
%         'm', 'LineWidth', 2);
% %     plot(radar_OIB.depth, median(radar_OIB.ages(:,OIB_near,:), 3)...
% %         + 2*std(squeeze(radar_OIB.ages(:,OIB_near,:)), [], 2), 'm--')
% %     plot(radar_OIB.depth, median(radar_OIB.ages(:,OIB_near,:), 3)...
% %         - 2*std(squeeze(radar_OIB.ages(:,OIB_near,:)), [], 2), 'm--')
%     h3 = plot(cores.(name).depth, mean(cores.(name).ages, 2), 'b', 'LineWidth', 2);
% %     plot(cores.(name).depth, mean(cores.(name).ages, 2) + ...
% %         2*std(cores.(name).ages, [], 2), 'b--')
% %     plot(cores.(name).depth, mean(cores.(name).ages, 2) - ...
% %         2*std(cores.(name).ages, [], 2), 'b--')
%     try
%         h4 = plot(radar_man.depth, radar_man.ages(:,man_near), 'k', 'LineWidth', 2);
%         legend([h1 h2 h3 h4], 'SEAT traces', 'OIB traces', 'Core', 'SEAT manual')
%     catch
%         legend([h1 h2 h3],'SEAT traces', 'OIB traces', 'Core')
%     end
%     title(strcat(text_name, ' age-depth scale'))
%     ylabel('Calendar years')
%     xlabel('Depth (m)')
%     hold off
%     
% %     % Save age-depth figure
% %     f1_name = strcat(name, '_age');
% %     export_fig(f1, fullfile(out_dir, f1_name), '-png');
%     
% %     % Plot of annual SMB estimates from OIB, SEAT core site, manual counts,
% %     % and core data (all realizations)
% %     f3 = figure('Position', [200 200 1300 700]);
% %     hold on
% %     title(strcat(text_name, ' annual SMB'))
% %     for n = OIB_idx
% %         h0 = plot(radar_OIB.SMB_yr{n}, median(radar_OIB.SMB{n}, 2), 'm', 'LineWidth', 0.5);
% %         h0.Color(4) = 0.02;
% %     end
% %     for n = 1:length(radar_i.SMB)
% %         h0 = plot(radar_i.SMB_yr{n}, median(radar_i.SMB{n}, 2), 'r', 'LineWidth', 0.5);
% %         h0.Color(4) = 0.02;
% %     end
% %     yr_end = min(cellfun(@length, radar_i.SMB_yr));
% %     SMB_data = cell2mat(cellfun(@(x) x(1:yr_end), radar_i.SMB, 'Uniform', false)')';
% %     h1 = plot(radar_i.SMB_yr{1}(1:yr_end), median(SMB_data, 2), 'r--', 'LineWidth', 1);
% %     yr_end = min(cellfun(@length, radar_OIB.SMB_yr(OIB_idx)));
% %     SMB_data = cell2mat(cellfun(@(x) x(1:yr_end), radar_OIB.SMB(OIB_idx), 'Uniform', false)')';
% %     h2 = plot(radar_OIB.SMB_yr{OIB_idx(1)}(1:yr_end), median(SMB_data, 2), 'm--', 'LineWidth', 1);
% %     
% %     h3 = plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2), 'b', 'LineWidth', 2);
% %     plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) + ...
% %         2*std(cores.(name).SMB, [], 2), 'b--');
% %     plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) - ...
% %         2*std(cores.(name).SMB, [], 2), 'b--');
% %     legend([h1 h2 h3], 'SEAT traces', 'OIB traces', 'Core')
% %     xlabel('Calendar years')
% %     ylabel('Annual SMB (mm w.e.)')
% %     f3_name = strcat(name, '_SMB');
% %     %     export_fig(f1, strcat(out_dir, f3_name), '-png');
% %     %     close(f3)
%     
%     % Plot of annual SMB estimates from OIB, SEAT core site, manual counts,
%     % and core data (only nearest trace)
%     f2 = figure('Position', [200 200 1000 700]);
%     hold on
%     title(strcat(text_name, ' annual SMB'))
%         for n = 1:Ndraw
%             h0 = plot(radar_OIB.SMB_yr{OIB_near}, radar_OIB.SMB{OIB_near}(:,n), 'm', 'LineWidth', 0.5);
%             h0.Color(4) = 0.05;
%         end
%         for n = 1:Ndraw
%             h0 = plot(radar_i.SMB_yr{SEATi_near}, radar_i.SMB{SEATi_near}(:,n), 'r', 'LineWidth', 0.5);
%             h0.Color(4) = 0.05;
%         end
%         for n = 1:Ndraw
%             h0 = plot(cores.(name).SMB_yr, cores.(name).SMB(:,n), 'b', 'LineWidth', 0.5);
%             h0.Color(4) = 0.05;
%         end
%     h1 = plot(radar_i.SMB_yr{SEATi_near}, median(radar_i.SMB{SEATi_near}, 2),...
%         'r', 'LineWidth', 2);
% %     plot(radar_i.SMB_yr{SEATi_near}, median(radar_i.SMB{SEATi_near}, 2) + ...
% %         2*std(radar_i.SMB{SEATi_near}, [], 2), 'r--');
% %     plot(radar_i.SMB_yr{SEATi_near}, median(radar_i.SMB{SEATi_near}, 2) - ...
% %         2*std(radar_i.SMB{SEATi_near}, [], 2), 'r--');
%     h2 = plot(radar_OIB.SMB_yr{OIB_near}, median(radar_OIB.SMB{OIB_near}, 2), ...
%         'm', 'LineWidth', 2);
% %     plot(radar_OIB.SMB_yr{OIB_near}, median(radar_OIB.SMB{OIB_near}, 2) + ...
% %         2*std(radar_OIB.SMB{OIB_near}, [], 2), 'm--');
% %     plot(radar_OIB.SMB_yr{OIB_near}, median(radar_OIB.SMB{OIB_near}, 2) - ...
% %         2*std(radar_OIB.SMB{OIB_near}, [], 2), 'm--');
%     h3 = plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2), ...
%         'b', 'LineWidth', 2);
% %     plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) + ...
% %         2*std(cores.(name).SMB, [], 2), 'b--');
% %     plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) - ...
% %         2*std(cores.(name).SMB, [], 2), 'b--');
%     xlim([min([min(radar_i.SMB_yr{SEATi_near}) min(radar_OIB.SMB_yr{OIB_near}) ...
%         min(cores.(name).SMB_yr)]) ...
%         max([max(radar_i.SMB_yr{SEATi_near}) max(radar_OIB.SMB_yr{OIB_near}) ...
%         max(cores.(name).SMB_yr)])])
%     legend([h1 h2 h3], 'SEAT traces', 'OIB traces', 'Core')
%     xlabel('Calendar years')
%     ylabel('Annual SMB (mm w.e.)')
%     hold off
%     
% %     % Save SMB figure
% %     f2_name = strcat(name, '_SMB');
% %     export_fig(f2, fullfile(out_dir, f2_name), '-png');
%     
% %     res = zeros(size(radar_i.Easting));
% %     for n = 1:length(radar_i.Easting)
% %         accum_n = mean(mean(accumulation_data(radar_i.Easting(n), radar_i.Northing(n))));
% %         res(n) = mean(mean(radar_i.SMB{n}, 2)) - accum_n;
% %     end
% %     
% %     res_core = mean(mean(cores.SEAT10_4.SMB, 2)) - accum_n;
% %     
% %     figure
% %     title(strcat(text_name, ' mean accum bias'))
% %     hold on
% %     h1 = histogram(res, 50);
% %     h2 = histogram(res_core);
% %     xlabel('mm/a bias (relative to Arthern et al 2006)')
% %     legend([h1 h2], 'SEAT radar', 'Core')
% %     hold off
%     
% 
%     %%% Welch's tests bewteen SEAT and OIB results
%     p_vals = zeros(min([length(radar_i.depth) length(radar_OIB.depth) ...
%         length(cores.(name).depth)]),1);
%     for k = 1:length(p_vals)
% %         [~, p_vals(k)] = ttest2(radar_i.ages(k,SEATi_near,:), radar_OIB.ages(k,OIB_near,:));
%         [~, p_vals(k)] = ttest2(cores.(name).ages(k,:), ...
%             squeeze(radar_OIB.ages(k,OIB_near,:)));
%     end
%     
%     
%     
% %     %%% trend tests
% %     yrs_endpts = [min([max(cores.(name).SMB_yr) max(radar_i.SMB_yr{SEATi_near}) ...
% %         max(radar_OIB.SMB_yr{OIB_near})]) max([min(cores.(name).SMB_yr) ...
% %         min(radar_i.SMB_yr{SEATi_near}) min(radar_OIB.SMB_yr{OIB_near})])];
% %     
% %     core_idxn = find(cores.(name).SMB_yr == yrs_endpts(1)):...
% %         find(cores.(name).SMB_yr == yrs_endpts(2));
% %     SEATi_idxn = find(radar_i.SMB_yr{SEATi_near} == yrs_endpts(1)):...
% %         find(radar_i.SMB_yr{SEATi_near} == yrs_endpts(2));
% %     OIB_idxn = find(radar_OIB.SMB_yr{OIB_near} == yrs_endpts(1)):...
% %         find(radar_OIB.SMB_yr{OIB_near} == yrs_endpts(2));
% %     
% %     yr_idx = [core_idxn' SEATi_idxn' OIB_idxn'];
%     
% %     trends.(name) = zeros(Ndraw, 3);
% %     p_vals.(name) = zeros(Ndraw, 3);
% %     for n = 1:Ndraw
% %         [b, stats] = robustfit(cores.(name).SMB_yr(yr_idx(:,1)), ...
% %             cores.(name).SMB(yr_idx(:,1),n));
% %         trends.(name)(n,1) = b(2);
% %         p_vals.(name)(n,1) = mean(stats.p);
% %         
% %         [b, stats] = robustfit(radar_i.SMB_yr{SEATi_near}(yr_idx(:,2)), ...
% %             radar_i.SMB{SEATi_near}(yr_idx(:,2),n));
% %         trends.(name)(n,2) = b(2);
% %         p_vals.(name)(n,2) = mean(stats.p);
% %         
% %         [b, stats] = robustfit(radar_OIB.SMB_yr{OIB_near}(yr_idx(:,3)), ...
% %             radar_OIB.SMB{OIB_near}(yr_idx(:,3),n));
% %         trends.(name)(n,3) = b(2);
% %         p_vals.(name)(n,3) = mean(stats.p);
% %     end
% %     
% %     sig = sum(p_vals.(name)<=0.05)/Ndraw;
% 
%     
% end
% 
% 
% 
% %% Comparison of SEAT and OIB estimates
% 
% dist_SEAT = pdist2([radar_SEAT.Easting' radar_SEAT.Northing'], ...
%     [radar_OIB.Easting' radar_OIB.Northing']);
% % [~, dist_idx] = min(dist_OIB, [], 2);
% [~, dist_idx] = min(dist_SEAT, [], 2);
% 
% res = zeros(size(radar_SEAT.data_smooth));
% for i = 1:length(radar_SEAT.Easting)
%     res(:,i) = median(squeeze(radar_SEAT.ages(:,i,:)), 2) - ...
%         median(squeeze(radar_OIB.ages(:,dist_idx(i),:)), 2) + 1;
% %     ho = plot(radar_SEAT.depth, res(:,i), 'r');
% %     ho.Color(4) = 0.01;
% end
% 
% % Find the root median squared error
% RMSE = sqrt(median(res.^2, 2));
% % MAE = median(abs(res), 2);
% p = polyfit(radar_SEAT.depth,RMSE, 1);
% bias = p(1);
% 
% figure
% histogram(res(end,:), 100)
% 
% % res_med = median(res, 2);
% % res_std = std(res, [], 2);
% % figure
% % hold on
% % plot(radar_SEAT.depth, res_med, 'k')
% % plot(radar_SEAT.depth, res_med + 2*res_std, 'k--')
% % plot(radar_SEAT.depth, res_med - 2*res_std, 'k--')
% % hold off
% 
% %%
% radar = radar_OIB;
% 
% % Calculate mean accumulation rate and std at each location
% SMB_mean = cellfun(@mean, radar.SMB, 'UniformOutput', 0);
% SMB_std = cellfun(@std, radar.SMB, 'UniformOutput', 0);
% 
% trend = cell(1, length(radar.SMB));
% p_vals = cell(1, length(radar.SMB));
% for i = 1:length(radar.SMB)
%     
%     trend_i = zeros(1, Ndraw);
%     p_val_i = zeros(1, Ndraw);
%     for j = 1:Ndraw
%         % Regression using regress function (average of MC simulations)
%         [b, ~, ~, ~, stats] = regress(radar.SMB{i}(:,j), ...
%             [ones(length(radar.SMB_yr{i}), 1) radar.SMB_yr{i}]);
%         trend_i(j) = b(2);
%         p_val_i(j) = stats(3);
%     end
%     trend{i} = trend_i;
%     p_vals{i} = p_val_i;
% end
% 
% trend_mean = cellfun(@mean, trend);
% trend_std = cellfun(@std, trend);
% 
% % Find indices of significant trends (at 95% confidence level)
% p_star = cellfun(@(x) find(x<=0.05), p_vals, 'UniformOutput', 0);
% 
% % Calculate percentage of MC realizations with significant trends
% p_ratio = cellfun(@(sig, all) length(sig)/length(all), p_star, p_vals);
% 
% % Regression of SMB trend on mean SMB (all traces)
% [trend_b, trend_int, ~, ~, trend_stats] = regress(trend_mean', ...
%     [ones(length(trend_mean), 1) cellfun(@mean, SMB_mean)']);
% [p_b, S_b] = polyfit(cellfun(@mean, SMB_mean), trend_mean, 1);
% [trend_Y, trend_D] = polyconf(p_b, cellfun(@mean, SMB_mean), S_b);
% 
% %%
% 
% subset = 1:25:length(SMB_mean);
% 
% % Boxplots of SMB mean
% mean_box = (cell2mat(SMB_mean'))';
% 
% f6 = figure('Position', [50 550 1800 450]);
% % boxplot([(mean(cores.(name).SMB))' mean_box(:,subset)])
% boxplot(mean_box(:,subset))
% ylabel('Mean SMB (mm/a)')
% set(gca,'Xticklabel',[])
% title('Mean SMB')
% 
% f6_name = 'SMB_box';
% % export_fig(f6, fullfile(out_dir, f6_name), '-png');
% % close(f6)
% 
% % Boxplots of trends in SMB
% trend_box = (cell2mat(trend'))';
% 
% f7 = figure('Position', [50 50 1800 450]);
% title('SMB trends')
% % boxplot([trend_core' trend_box(:,subset)])
% boxplot(trend_box(:,subset))
% ylabel('SMB trend (mm/a)')
% set(gca,'Xticklabel',[])
% title('SMB trends')
% 
% f7_name = strcat(name, '_trend_box');
% % export_fig(f7, fullfile(out_dir, f7_name), '-png');
% % close(f7)
% 
% %% Compare mean accumulation to Arthern et al 2006
% %%% Can also add in White et al comparison using Favier compilation
% 
% [accum_E, accum_N, accum_A] = accumulation_data(radar_OIB.Easting, ...
%     radar_OIB.Northing, 'xy');
% accum_OIB = cellfun(@(x) median(median(x,2)), radar_OIB.SMB);
% 
% figure
% hold on
% imagesc(accum_E(1,:), accum_N(:,1), accum_A)
% scatter(radar_OIB.Easting, radar_OIB.Northing, 100, accum_OIB, 'filled')
% colorbar
% hold off
% 
% figure
% hold on
% plot(radar_OIB.Easting, accum_OIB, 'm')
% plot(accum_E(2,:), accum(2,:), 'c')
% 
% % res = zeros(size(radar_OIB.Easting));
% % for i = 1:length(radar_OIB.Easting)
% %     accum_i = mean(mean(accumulation_data(radar_OIB.Easting(i), radar_OIB.Northing(i))));
% %     res(i) = median(median(radar_OIB.SMB{i}, 2)) - accum_i;
% % end
% % % res_core = mean(mean(cores.SEAT10_4.SMB, 2)) - accum_i;
% % figure
% % plot(radar_OIB.Easting, res, 'r.')

