


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
% SEAT_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
%     'SEAT2010Kuband/SEAT10_4toSEAT10_6/SMB_results/', wild));
SEAT_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
    'SEAT2010Kuband/allSEAT10_4toSEAT10_6/SMB_results/', wild));

seat_E = [];
seat_N = [];
seat_SMB_MC = [];
seat_yr = [];
for i = 1:length(SEAT_files)
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'Easting');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'Northing');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'SMB');
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'SMB_yr');
    
    seat_E = [seat_E Easting];
    seat_N = [seat_N Northing];
    seat_SMB_MC = [seat_SMB_MC SMB];
    seat_yr = [seat_yr SMB_yr];
end
seat_SMB = cellfun(@(x) mean(x, 2), seat_SMB_MC, 'UniformOutput', 0);
seat_std = cellfun(@(x) std(x, [], 2), seat_SMB_MC, 'UniformOutput', 0);


wild = '*.mat';
OIB_files = dir(fullfile(data_path, ...
    'IceBridge/SEAT10_4to10_6/2011_SNO/SMB_results', wild));

oib_E = [];
oib_N = [];
oib_elev = [];
oib_SMB_MC = [];
oib_yr = [];
for i = 1:length(OIB_files)
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Easting');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Northing');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB_yr');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'elev');
    
    oib_E = [oib_E Easting];
    oib_N = [oib_N Northing];
    oib_elev = [oib_elev, elev];
    oib_SMB_MC = [oib_SMB_MC SMB];
    oib_yr = [oib_yr SMB_yr];
end
oib_SMB = cellfun(@(x) mean(x, 2), oib_SMB_MC, 'UniformOutput', 0);
oib_std = cellfun(@(x) std(x, [], 2), oib_SMB_MC, 'UniformOutput', 0);
clear Easting Northing SMB SMB_yr i wild

%%

yr_start = 1980;
yr_end = 2009;
year = (yr_end:-1:yr_start)';

seat_idx = cellfun(@(x) max(x)>=yr_end && min(x)<=yr_start, seat_yr);
SEAT_E = seat_E(seat_idx);
SEAT_N = seat_N(seat_idx);
SEAT_SMB = seat_SMB(seat_idx);
SEAT_yr = seat_yr(seat_idx);
oib_idx = cellfun(@(x) max(x)>=yr_end && min(x)<=yr_start, oib_yr);
OIB_E = oib_E(oib_idx);
OIB_N = oib_N(oib_idx);
OIB_elev = oib_elev(oib_idx);
OIB_SMB = oib_SMB(oib_idx);
OIB_yr = oib_yr(oib_idx);

SEAT_start = cellfun(@(x) find(x==yr_start, 1), SEAT_yr, 'UniformOutput', 0);
SEAT_end = cellfun(@(x) find(x==yr_end, 1), SEAT_yr, 'UniformOutput', 0);
mSEAT_SMB = cellfun(@(x,y,z) x(y:z), SEAT_SMB, SEAT_end, SEAT_start, ...
    'UniformOutput', 0);
% [coeff,~,~,~,stats] = cellfun(@(x) regress(x, [ones(length(year),1) year]),...
%     mSEAT_SMB, 'UniformOutput', 0);
% SEAT_beta = cellfun(@(x) x(2), coeff);
% SEAT_pval = cellfun(@(x) x(3), stats);
[coeff, stats] = cellfun(@(x) robustfit(year, x), mSEAT_SMB, 'UniformOutput', 0);
SEAT_beta = cellfun(@(x) x(2), coeff);
SEAT_pval = cellfun(@(x) x.p(2), stats);

OIB_start = cellfun(@(x) find(x==yr_start, 1), OIB_yr, 'UniformOutput', 0);
OIB_end = cellfun(@(x) find(x==yr_end, 1), OIB_yr, 'UniformOutput', 0);
mOIB_SMB = cellfun(@(x,y,z) x(y:z), OIB_SMB, OIB_end, OIB_start, ...
    'UniformOutput', 0);
% [coeff,~,~,~,stats] = cellfun(@(x) regress(x, [ones(length(year),1) year]),...
%     mOIB_SMB, 'UniformOutput', 0);
% OIB_beta = cellfun(@(x) x(2), coeff);
% OIB_pval = cellfun(@(x) x(3), stats);
[coeff, stats] = cellfun(@(x) robustfit(year, x), mOIB_SMB, 'UniformOutput', 0);
OIB_beta = cellfun(@(x) x(2), coeff);
OIB_pval = cellfun(@(x) x.p(2), stats);


my_cores = {'SEAT10_4', 'SEAT10_5', 'SEAT10_6'};
cores_E = zeros(1, length(my_cores));
cores_N = zeros(1, length(my_cores));
cores_SMB = zeros(length(year), length(my_cores));
cores_beta = zeros(1, length(my_cores));
cores_pval = zeros(1, length(my_cores));
for k = 1:length(my_cores)
    core_k = cores.(my_cores{k});
    cores_E(k) = core_k.Easting;
    cores_N(k) = core_k.Northing;
    core_start = find(core_k.SMB_yr==yr_end, 1);
    core_end = find(core_k.SMB_yr==yr_start, 1);
    SMB_mean = mean(core_k.SMB, 2);
    SMB_k = SMB_mean(core_start:core_end);
    cores_SMB(:,k) = SMB_k;
%     [coeff,~,~,~,stats] = regress(SMB_k, [ones(length(year),1) year]);
%     cores_beta(k) = coeff(2);
%     cores_pval(k) = stats(3);
    [coeff, stats] = robustfit(year, SMB_k);
    cores_beta(k) = coeff(2);
    cores_pval(k) = stats.p(2);
end

figure
hold on
plot(SEAT_E, mean(cell2mat(mSEAT_SMB)), 'r.')
plot(OIB_E, mean(cell2mat(mOIB_SMB)), 'm.')
plot(cores_E, mean(cores_SMB), 'bx', 'MarkerSize', 10, 'LineWidth', 3)
xlabel('Position (Easting)')
ylabel(['Mean SMB (' num2str(yr_start) '-' num2str(yr_end) ')'])
legend('Location', 'northeast', 'SEAT', 'OIB', 'cores')
hold off

subset = 1:20:length(mOIB_SMB);
figure
boxplot(fliplr(cell2mat(mOIB_SMB(subset))), 'PlotStyle', 'compact')
ylabel(['OIB annual SMB' num2str(yr_start) '-' num2str(yr_end) ' (mm/a)'])
a = get(get(gca, 'children'), 'children');
set(a, 'Color', 'm')
% hold on
% boxplot(cell2mat(mSEAT_SMB(1:20:length(mSEAT_SMB))), 'PlotStyle', 'compact')
set(gca, 'xtick', [])
set(gca, 'xticklabel', [])
xlabel('Easting')

OIB_sig = OIB_pval<=0.05;
cores_sig = cores_pval<=0.05;
figure
hold on
plot(OIB_E, OIB_beta, 'k.', 'MarkerSize', 5)
plot(OIB_E(OIB_sig), OIB_beta(OIB_sig), 'm.', 'MarkerSize', 5)
plot(cores_E, cores_beta, 'kx', 'MarkerSize', 15, 'LineWidth', 3)
plot(cores_E(cores_sig), cores_beta(cores_sig), 'bx', 'MarkerSize', 15, ...
    'LineWidth', 3)
xlabel('Trace position (Easting)')
ylabel(['Annual SMB trend' num2str(yr_start) '-' num2str(yr_end) ' (mm/a)'])
legend('Location', 'southeast', 'OIB (all)', 'OIB (significant)', ...
    'cores (all)', 'cores (sigificant)')
hold off

SEAT_sig = SEAT_pval<=0.05;
figure
hold on
plot(SEAT_E, SEAT_beta, 'ko')
plot(SEAT_E(SEAT_sig), SEAT_beta(SEAT_sig), 'ro')
plot(cores_E, cores_beta, 'kx', 'MarkerSize', 15, 'LineWidth', 3)
plot(cores_E(cores_sig), cores_beta(cores_sig), 'bx', 'MarkerSize', 15, ...
    'LineWidth', 3)
xlabel('Trace position (Easting)')
ylabel(['Annual SMB trend' num2str(yr_start) '-' num2str(yr_end) ' (mm/a)'])
legend('Location', 'southeast', 'SEAT (all)', 'SEAT (significant)', ...
    'cores (all)', 'cores (sigificant)')
hold off


% figure
% hold on
% % plot(SEAT_E, SEAT_beta./mean(cell2mat(mSEAT_SMB)), 'ko')
% % plot(SEAT_E(SEAT_sig), SEAT_beta(SEAT_sig)./...
% %     mean(cell2mat(mSEAT_SMB(SEAT_sig))), 'ro')
% plot(OIB_E, OIB_beta./mean(cell2mat(mOIB_SMB)), 'k.')
% plot(OIB_E(OIB_sig), OIB_beta(OIB_sig)./...
%     mean(cell2mat(mOIB_SMB(OIB_sig))), 'm.')
% plot(cores_E, cores_beta./mean(cores_SMB), 'kx', ...
%     'MarkerSize', 15, 'LineWidth', 3)
% plot(cores_E(cores_sig), cores_beta(cores_sig)./mean(cores_SMB(:,cores_sig)),...
%     'bx', 'MarkerSize', 15, 'LineWidth', 3)
% hold off


%% Compare mean accumulation to Arthern et al 2006
%%% Can also add in White et al comparison using Favier compilation

labels = strrep(cores.name, '_', '-');
basins = shaperead(strcat(data_path, ...
    'DEMs/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'));
Easting_lims = [min(cores.Easting)-10000 max(cores.Easting)+5000];
Northing_lims = [min(cores.Northing)-5000 max(cores.Northing)+5000];


[Arth_E, Arth_N, Arth_accum] = accumulation_data(Easting_lims, Northing_lims, 'xy');


figure('Position', [10 10 1400 800])
h0 = image(Arth_E(1,:), (Arth_N(:,1))', Arth_accum, 'CDataMapping', 'scaled');
set(gca, 'Ydir', 'normal')
hold on
h1 = plot(cores.Easting(1), cores.Northing(1), 'k', 'LineWidth', 2);
h2 = mapshow(basins, 'FaceAlpha', 0, 'LineWidth', 3);
h3 = scatter(OIB_E, OIB_N, 25, mean(cell2mat(mOIB_SMB)));
h4 = scatter(cores.Easting, cores.Northing, 50, 'b', 'filled');
text(cores.Easting, cores.Northing, strcat('\leftarrow', labels), ...
    'FontSize', 15, 'Interpreter', 'tex');
c0 = colorbar;
c0.Label.String = ['Mean annual SMB ' num2str(yr_start) '-' num2str(yr_end) ' (mm/a)'];
c0.Label.FontSize = 18;
graticuleps(-81:0.5:-77,-125:2:-105, 'c')
xlim(Easting_lims)
ylim(Northing_lims)
scalebarps
box on
mapzoomps('ne', 'insetsize', 0.30)
% legend([h1 h3 h4], 'WAIS Divide', 'OIB', 'firn cores',...
%     'Location', 'northwest')
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 18)
hold off


%% Elevation investigations

elev = cryosat2_interp(OIB_E, OIB_N);

figure
hold on
plot(OIB_E, elev, 'm')
plot(cores_E, cores.elev(3:5), 'bo')

figure
plot(elev, OIB_beta, 'm.')

%%

% % Calculate mean accumulation rate and std at each location
% SMB_med = cellfun(@median, radar.SMB, 'UniformOutput', 0);
% SMB_std = cellfun(@std, radar.SMB, 'UniformOutput', 0);
% 
% trend = cell(1, length(radar.SMB));
% p_val = cell(1, length(radar.SMB));
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
%     p_val{i} = p_val_i;
% end
% 
% trend_med = cellfun(@median, trend);
% trend_std = cellfun(@std, trend);
% 
% % Find indices of significant trends (at 95% confidence level)
% p_star = cellfun(@(x) find(x<=0.05), p_val, 'UniformOutput', 0);
% 
% % Calculate percentage of MC realizations with significant trends
% p_ratio = cellfun(@(sig, all) length(sig)/length(all), p_star, p_val);
% 
% % Regression of SMB trend on mean SMB (all traces)
% [trend_b, trend_int, ~, ~, trend_stats] = regress(trend_med', ...
%     [ones(length(trend_med), 1) cellfun(@mean, SMB_med)']);
% [p_b, S_b] = polyfit(cellfun(@mean, SMB_med), trend_med, 1);
% [trend_Y, trend_D] = polyconf(p_b, cellfun(@mean, SMB_med), S_b);



%% Diagnostic plots for bulk radar file

% % Compare significance of differences between core SMB and radar SMB
% yr_start = min([core_near1.SMB_yr(1) cellfun(@(x) x(1), radar.SMB_yr)]);
% yr_end = max([core_near1.SMB_yr(end) cellfun(@(x) x(end), radar.SMB_yr)]);
% yr_endpts = [yr_start yr_end];
% Ryr_idx = [find(radar.SMB_yr{i}==yr_endpts(1)) ...
%     find(radar.SMB_yr{i}==yr_endpts(2))];
% Cyr_idx = [find(core_near1.SMB_yr==yr_start) find(core_near1.SMB_yr==yr_end)];
% rSMB = cellfun(@(x) x(Ryr_idx(1):Ryr_idx(2),:), ...
%     radar.SMB, 'UniformOutput', 0);
% cSMB = core_near1.SMB(Cyr_idx(1):Cyr_idx(2),:);
% 
% rSMB_mu = cell2mat(cellfun(@(x) mean(x,2), rSMB, 'UniformOutput', 0));
% rSMB_SE = cell2mat(cellfun(@(x) std(x,[],2)/sqrt(Ndraw), rSMB, 'UniformOutput', 0));
% cSMB_mu = mean(cSMB,2);
% cSMB_SE = std(cSMB,[],2)/sqrt(Ndraw);
% SMB_CIup = (rSMB_mu-cSMB_mu) + 1.96*sqrt(rSMB_SE.^2 + cSMB_SE.^2);
% SMB_CIlow = (rSMB_mu-cSMB_mu) - 1.96*sqrt(rSMB_SE.^2 + cSMB_SE.^2);
% SMB_sig = SMB_CIup.*SMB_CIlow >=0;


% yr_start = min([core_near1.SMB_yr(1) cellfun(@(x) x(1), radar.SMB_yr)]);
% yr_end = max([core_near1.SMB_yr(end) cellfun(@(x) x(end), radar.SMB_yr)]);
% yr_endpts = [yr_start yr_end];
% Rendpts_idx = [find(radar.SMB_yr{i}==yr_endpts(1)) ...
%     find(radar.SMB_yr{i}==yr_endpts(2))];
% Cendpts_idx = [find(core_near1.SMB_yr==yr_start) find(core_near1.SMB_yr==yr_end)];
% rSMB = cellfun(@(x) x(Rendpts_idx(1):Rendpts_idx(2),:), ...
%     radar.SMB, 'UniformOutput', 0);
% SMB_res = cell2mat(cellfun(@(x) x - core_near1.SMB(Cendpts_idx(1):Cendpts_idx(2),:),...
%     rSMB, 'UniformOutput', 0));
% SMB_med = median(SMB_res,2);
% SMB_CI = zeros(size(SMB_res,1), 2);
% SMB_CI(:,1) = quantile(SMB_res, 0.95, 2);
% SMB_CI(:,2) = quantile(SMB_res, 0.05, 2);
% figure
% hold on
% plot(yr_start:-1:yr_end, SMB_med, 'r', 'LineWidth', 2)
% plot(yr_start:-1:yr_end, SMB_CI, 'b', 'LineWidth', 2)
% rline = refline(0,0);
% rline.LineStyle = ':';
% rline.Color = 'k';
% rline.LineWidth = 3;
% xlabel('Calendar year')
% ylabel('SMB bias')
% legend('Median bias', '0.95 CI')
% hold off
% % figure('Position', [200 200 1200 500])
% % boxplot(SMB_res', yr_start:-1:yr_end)
% % hold on
% % rline = refline(0,0);
% % rline.LineStyle = ':';
% % rline.Color = 'k';
% % rline.LineWidth = 2;
% % xlabel('Calendar year')
% % ylabel('SMB bias')
% % hold off





