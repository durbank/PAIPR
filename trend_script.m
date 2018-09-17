


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

% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Number of Monte Carlo simulations
Ndraw = 100;

% Import firn core data
% [cores] = import_cores(fullfile(data_path, 'Ice-cores/SEAT_cores/', ...
%     'DGK_core_data.xlsx'), Ndraw);

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

%% Define radar files to import/process

% radar_dir = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband/', ...
%     'ProcessedSEAT2010/transectSEAT10_4_5/');
% 
% % List all files matching 'wild' within radar directory
% wild = '*.mat';
% files = dir(fullfile(radar_dir, wild));
% 
% i = randi(length(files));
% file = strcat(radar_dir, files(i).name);

% Path to full SEAT transect
file = fullfile(data_path, 'radar/SEAT_Traverses/core-site_tests/', ...
    'layers_ku_band_SEAT10_5.mat');

% file = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband/', ...
%     'layers_ku_band_gridSEAT10_4.mat');

% % Path of the OIB file to process
% % SEAT10_4
% file = strcat(data_path, 'IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_272.nc');
% file = strcat(data_path, 'IceBridge/Snow Radar/2016/IRSNO1B_20161109_02_381.nc');
% file = strcat(data_path, 'IceBridge/Kuband/2016/IRKUB1B_20161109_02_381.nc');
% % SEAT10_5
% file = strcat(data_path, 'IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_257.nc');
% % SEAT10_6
% file = strcat(data_path, 'IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_242.nc');

%%

% for i = 1:length(files)

% Calculate radar ages and associated other data
% [radar] = radar_age(file, cores, Ndraw);
[radar] = radar_RT(file, cores, Ndraw);

% Calculate annual accumulation rates from data
[radar] = calc_SWE(radar, Ndraw);


% clip = round (3000/25);
% radar0 = radar;
% radar = struct('collect_date', radar.collect_date, 'Easting', radar.Easting(clip:end-clip),...
%     'Northing', radar.Northing(clip:end-clip), 'dist', radar.dist(clip:end-clip),...
%     'depth', radar.depth, 'data_smooth', radar.data_smooth(:,clip:end-clip),...
%     'peaks', radar.peaks(:,clip:end-clip), 'groups', radar.groups(:,clip:end-clip),...
%     'likelihood', radar.likelihood(:,clip:end-clip), 'ages', radar.ages(:,clip:end-clip,:));
% radar.SMB_yr =  radar0.SMB_yr(clip:end-clip);
% radar.SMB = radar0.SMB(clip:end-clip);

% output_path = fullfile(data_path, 'radar/SEAT_Traverses/results_data/gridSEAT10_4.mat');
% save(output_path, '-struct', 'radar', '-v7.3')

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


%% Diagnostic plots for single random trace

% Trace idx to investigate
i = randi(length(radar.SMB));

% Find the nearest cores to the radar data (for comparison plots)
[~, cores_near_idx] = sort(pdist2([radar.Easting(i) radar.Northing(i)], ...
    [cores.Easting' cores.Northing'], 'Euclidean'));
core_near1 = cores.(cores.name{cores_near_idx(1)});
core_near2 = cores.(cores.name{cores_near_idx(2)});
core_near3 = cores.(cores.name{cores_near_idx(3)});

% Calculate the median age-depth scale and std for radar trace i
age_med = median(squeeze(radar.ages(:,i,:)), 2);
age_std = std(squeeze(radar.ages(:,i,:)), [], 2);

% Plot full radargram
yr_idx = logical([diff(floor(age_med)); 0]);
layer_idx = radar.likelihood(:,i) >= 0.10;
% depth = radar.depth(yr_idx);
depth = radar.depth(layer_idx);
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
xlim([radar.dist(1) radar.dist(end)])
ylim([0 radar.depth(end)])
set(gca, 'Ydir', 'reverse', 'FontSize', 18)
hold off

% Age-depth scale comparison between radar trace and nearest cores
figure
hold on
h1 = plot(core_near1.depth, median(core_near1.ages, 2), 'b', 'LineWidth', 2);
plot(core_near1.depth, median(core_near1.ages, 2) + std(core_near1.ages, [], 2), 'b--')
plot(core_near1.depth, median(core_near1.ages, 2) - std(core_near1.ages, [], 2), 'b--')
h2 = plot(core_near2.depth, core_near2.age, 'c', 'LineWidth', 2);
h3 = plot(core_near3.depth, core_near3.age, 'c--', 'LineWidth', 1);
h4 = plot(radar.depth, age_med, 'r', 'LineWidth', 2);
plot(radar.depth, age_med + age_std, 'r--', 'LineWidth', 0.5)
plot(radar.depth, age_med - age_std, 'r--', 'LineWidth', 0.5)
ylabel('Calendar Year')
xlabel('Depth (m)')
legend([h1 h2 h3 h4], 'Nearest core age (manual)', '2nd nearest core', ...
    '3rd nearest core', 'Radar age (automated)', 'Location', 'ne')
set(gca, 'FontSize', 10)
hold off

% Compare annual accumulation between ith radar trace and nearest 2 cores
% (along with estimated uncertainties for each)
figure
hold on
h1 = plot(core_near1.SMB_yr, median(core_near1.SMB, 2), 'b', 'LineWidth', 2);
plot(core_near1.SMB_yr, median(core_near1.SMB, 2) + std(core_near1.SMB, [], 2), 'b--')
plot(core_near1.SMB_yr, median(core_near1.SMB, 2) - std(core_near1.SMB, [], 2), 'b--')

h2 = plot(core_near2.SMB_yr, median(core_near2.SMB, 2), 'c', 'LineWidth', 2);
plot(core_near2.SMB_yr, median(core_near2.SMB, 2) + std(core_near2.SMB, [], 2), 'c--')
plot(core_near2.SMB_yr, median(core_near2.SMB, 2) - std(core_near2.SMB, [], 2), 'c--')

h3 = plot(radar.SMB_yr{i}, median(radar.SMB{i}, 2), 'r', 'LineWidth', 2);
plot(radar.SMB_yr{i}, median(radar.SMB{i}, 2) + std(radar.SMB{i}, [], 2), 'r--')
plot(radar.SMB_yr{i}, median(radar.SMB{i}, 2) - std(radar.SMB{i}, [], 2), 'r--')
legend([h1 h2 h3], 'Nearest firn core', '2nd nearest core', 'Ku radar')
xlabel('Calendar Year')
ylabel('Annual accumulation (mm w.e.)')
hold off

%% Diagnostic plots for bulk radar file

% Compare significance of differences between core ages and radar ages
ages_idx = [1 min([length(core_near1.depth) length(radar.depth)])];
rAGES = radar.ages(1:ages_idx(2),:,:);
% cAGES = repmat(reshape(core_near1.ages, size(core_near1.ages,1), 1,...
%     size(core_near1.ages,2)), [1 length(radar.dist), 1]);

rAGES_mu = mean(rAGES,3);
rAGES_SE = std(rAGES, [], 3)/sqrt(Ndraw);
cAGES_mu = mean(core_near1.ages, 2);
cAGES_SE = std(core_near1.ages,[],2)/sqrt(Ndraw);
ages_CIup = (rAGES_mu-cAGES_mu) + 1.96*sqrt(rAGES_SE.^2 + cAGES_SE.^2);
ages_CIlow = (rAGES_mu-cAGES_mu) - 1.96*sqrt(rAGES_SE.^2 + cAGES_SE.^2);
ages_sig = ages_CIup.*ages_CIlow >=0;


% ages_res = reshape(rAGES - cAGES, Rendpts_idx(2), size(cAGES,2)*size(cAGES,3));
% ages_med = median(ages_res,2);
% ages_CI = zeros(Rendpts_idx(2), 2);
% ages_CI(:,1) = quantile(ages_res, 0.95, 2);
% ages_CI(:,2) = quantile(ages_res, 0.05, 2);
% figure
% hold on
% plot(radar.depth(1:Rendpts_idx(2)), ages_med, 'r', 'LineWidth', 2)
% plot(radar.depth(1:Rendpts_idx(2)), ages_CI, 'b', 'LineWidth', 2)
% rline = refline(0,0);
% rline.LineStyle = ':';
% rline.Color = 'k';
% rline.LineWidth = 3;
% xlabel('Depth (m)')
% ylabel('Age bias')
% legend('Median bias', '0.95 CI')
% hold off
% % figure('Position', [200 200 1200 500])
% % boxplot(ages_res(1:50:end,1:75:end)', 0:0.02*50:radar.depth(Rendpts_idx(2)))
% % hold on
% % rline = refline(0,0);
% % rline.LineStyle = ':';
% % rline.Color = 'k';
% % rline.LineWidth = 2;
% % xlabel('Depth (m)')
% % ylabel('Age bias')
% % hold off


% Compare significance of differences between core SMB and radar SMB
yr_start = min([core_near1.SMB_yr(1) cellfun(@(x) x(1), radar.SMB_yr)]);
yr_end = max([core_near1.SMB_yr(end) cellfun(@(x) x(end), radar.SMB_yr)]);
yr_endpts = [yr_start yr_end];
Ryr_idx = [find(radar.SMB_yr{i}==yr_endpts(1)) ...
    find(radar.SMB_yr{i}==yr_endpts(2))];
Cyr_idx = [find(core_near1.SMB_yr==yr_start) find(core_near1.SMB_yr==yr_end)];
rSMB = cellfun(@(x) x(Ryr_idx(1):Ryr_idx(2),:), ...
    radar.SMB, 'UniformOutput', 0);
cSMB = core_near1.SMB(Cyr_idx(1):Cyr_idx(2),:);

rSMB_mu = cell2mat(cellfun(@(x) mean(x,2), rSMB, 'UniformOutput', 0));
rSMB_SE = cell2mat(cellfun(@(x) std(x,[],2)/sqrt(Ndraw), rSMB, 'UniformOutput', 0));
cSMB_mu = mean(cSMB,2);
cSMB_SE = std(cSMB,[],2)/sqrt(Ndraw);
SMB_CIup = (rSMB_mu-cSMB_mu) + 1.96*sqrt(rSMB_SE.^2 + cSMB_SE.^2);
SMB_CIlow = (rSMB_mu-cSMB_mu) - 1.96*sqrt(rSMB_SE.^2 + cSMB_SE.^2);
SMB_sig = SMB_CIup.*SMB_CIlow >=0;


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


% Mean SMB across entire radargram (mean of all realizations for each trace)
% figure
% scatter(radar.Easting, radar.Northing, 30, cellfun(@mean, SMB_mean), 'filled')
% hcb = colorbar;
% ylabel(hcb, 'Mean annual SMB (mm/a')
% colormap(cool)
% 
% % Mean trend in SMB across entire radargram
% figure
% scatter(radar.Easting, radar.Northing, 30, ...
%     100*trend_mean./cellfun(@mean, SMB_mean), 'filled')
% hcb = colorbar;
% ylabel(hcb, 'Trend in SMB (% of mean per year)')
% colormap(cool)


SMB_med = cellfun(@(x) median(x,2), radar.SMB, 'UniformOutput', 0);
[regress_coeff, stats] = cellfun(@robustfit, radar.SMB_yr, SMB_med, 'UniformOutput', 0);
SMB_trend = cellfun(@(x) x(2), regress_coeff);
trend_SE = cellfun(@(x) x.se(2), stats);

SMB_mean = cellfun(@mean, SMB_med);

% Mean SMB vs mean trend
figure
hold on
plot(SMB_mean, 100*SMB_trend./SMB_mean, 'b.')
plot([min(SMB_mean) max(SMB_mean)], ...
    [median(100*SMB_trend./SMB_mean) median(100*SMB_trend./SMB_mean)], 'r--')
xlabel('Mean annual accumulation (mm w.e.)')
ylabel('Accumulation trend (% of mean per year)')
hold off