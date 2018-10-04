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
% SEAT_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
%     'SEAT2010Kuband/SEAT10_4toSEAT10_6/SMB_results/', wild));
SEAT_files = dir(fullfile(data_path, 'radar/SEAT_Traverses/',...
    'SEAT2010Kuband/allSEAT10_4toSEAT10_6/SMB_results/', wild));

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
SEAT_SMB = cellfun(@(x) mean(x, 2), SEAT_SMB_MC, 'UniformOutput', 0);
SEAT_std = cellfun(@(x) std(x, [], 2), SEAT_SMB_MC, 'UniformOutput', 0);


wild = '*.mat';
OIB_files = dir(fullfile(data_path, ...
    'IceBridge/SEAT10_4to10_6/2011_SNO/SMB_results', wild));

OIB_E = [];
OIB_N = [];
OIB_SMB_MC = [];
OIB_yr = [];
for i = 1:length(OIB_files)
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Easting');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'Northing');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB');
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'SMB_yr');
    
    OIB_E = [OIB_E Easting];
    OIB_N = [OIB_N Northing];
    OIB_SMB_MC = [OIB_SMB_MC SMB];
    OIB_yr = [OIB_yr SMB_yr];
end
clear Easting Northing SMB SMB_yr i wild
OIB_SMB = cellfun(@(x) mean(x, 2), OIB_SMB_MC, 'UniformOutput', 0);
OIB_std = cellfun(@(x) std(x, [], 2), OIB_SMB_MC, 'UniformOutput', 0);



% SEAT_full = load(fullfile(SEAT_files(2).folder, SEAT_files(2).name));
% SEAT_layers = cell(1, max(SEAT_full.groups(:)));
% for k = 1:length(SEAT_layers)
%     SEAT_layers{k} = find(SEAT_full.groups(:,141:end)==k);
% end
% SEAT_layers = SEAT_layers(~cellfun(@isempty, SEAT_layers));
% 
% f1 = figure('Position', [50 50 1800 900]);
% imagesc(SEAT_full.Easting(141:end), SEAT_full.depth, ...
%     SEAT_full.data_smooth(:,141:end), [-2 2])
% hold on
% for k = 1:length(SEAT_layers)
%     [r,c] = ind2sub(size(SEAT_full.groups(:,141:end)), SEAT_layers{k});
%     plot(SEAT_full.Easting(c+141-1),SEAT_full.depth(r), 'r', 'LineWidth',2)
% end
% xlabel('Easting coordinate')
% ylabel('Depth (m)')
% hold off
% % % Save SEAT layers figure
% % f1_name = strcat('SEAT_layers');
% % out_dir = fullfile('F:/results/');
% % export_fig(f1, fullfile(out_dir, f1_name), '-png');
% 
% 
% OIB_full = load(fullfile(OIB_files(7).folder, OIB_files(7).name));
% OIB_layers = cell(1, max(OIB_full.groups(:)));
% for k = 1:length(OIB_layers)
%     OIB_layers{k} = find(OIB_full.groups(:,150:end)==k);
% end
% OIB_layers = OIB_layers(~cellfun(@isempty, OIB_layers));
% 
% f2 = figure('Position', [50 50 1800 900]);
% imagesc(OIB_full.Easting(150:end), OIB_full.depth, ...
%     OIB_full.data_smooth(:,150:end), [-2 2])
% hold on
% for k = 1:length(OIB_layers)
%     [r,c] = ind2sub(size(OIB_full.groups(:,150:end)), OIB_layers{k});
%     plot(OIB_full.Easting(c+150-1), OIB_full.depth(r), ...
%         'r', 'LineWidth',2)
% end
% xlabel('Easting coordinate')
% ylabel('Depth (m)')
% hold off
% % % Save OIB layers figure
% % f2_name = strcat('OIB_layers');
% % out_dir = fullfile('F:/results/');
% % export_fig(f2, fullfile(out_dir, f2_name), '-png');


%% Compare median SMB estimates within 10 km of site 4 (core site SEAT10-4)

cores_loop = {'SEAT10_4' 'SEAT10_5' 'SEAT10_6'};
near_dist = 6500;
SMB_sites = struct();
out_dir = fullfile('F:/results/');

for i = 1:length(cores_loop)
    
    core_i = cores.(cores_loop{i});
    
    D_SEATi = pdist2([core_i.Easting, core_i.Northing], [SEAT_E' SEAT_N']);
    SEATi_idx = D_SEATi <= near_dist;
    [~, SEATi_near] = min(D_SEATi);
    
    D_OIBi = pdist2([core_i.Easting, core_i.Northing], [OIB_E' OIB_N']);
    OIBi_idx = D_OIBi <= near_dist;
    [~, OIBi_near] = min(D_OIBi);
    
    yr_start = min([core_i.SMB_yr(1) cellfun(@(x) x(1), SEAT_yr(SEATi_idx)) ...
        cellfun(@(x) x(1), OIB_yr(OIBi_idx))]);
    yr_end = max([core_i.SMB_yr(end) cellfun(@(x) x(end), SEAT_yr(SEATi_idx)) ...
        cellfun(@(x) x(end), OIB_yr(OIBi_idx))]);
    siteI_yr = (yr_start:-1:yr_end)';
    
    SEAT_start = cellfun(@(x) find(x==yr_start, 1), SEAT_yr(SEATi_idx));
    SEAT_end = cellfun(@(x) find(x==yr_end, 1), SEAT_yr(SEATi_idx));
    OIB_start = cellfun(@(x) find(x==yr_start, 1), OIB_yr(OIBi_idx));
    OIB_end = cellfun(@(x) find(x==yr_end, 1), OIB_yr(OIBi_idx));
    
    SMB_tmp = SEAT_SMB(SEATi_idx);
    SEATi_SMB = zeros(length(siteI_yr), length(SEAT_start));
    for j = 1:length(SEAT_start)
        
        SEATi_SMB(:,j) = SMB_tmp{j}(SEAT_start(j):SEAT_end(j));
    end
    % SEAT4_SMB = movmean(SEAT4_SMB, 3);
    
    SMB_tmp = OIB_SMB(OIBi_idx);
    OIBi_SMB = zeros(length(siteI_yr), length(OIB_start));
    for j = 1:length(OIB_start)
        
        OIBi_SMB(:,j) = SMB_tmp{j}(OIB_start(j):OIB_end(j));
    end
    % OIB4_SMB = movmean(OIB4_SMB, 3);
    
    
    f1 = figure;
    hold on
    SEAT_SMB_near = SEAT_SMB_MC{SEATi_near}(SEAT_start(1):...
        SEAT_start(1)+length(siteI_yr)-1,:);
    for n = 1:size(SEAT_SMB_near, 2)
        h0 = plot(siteI_yr, SEAT_SMB_near(:,n), 'r', 'LineWidth', 0.5);
        h0.Color(4) = 0.02;
    end
    h1 = plot(siteI_yr, mean(SEAT_SMB_near, 2), 'r', 'LineWidth', 2);
    plot(siteI_yr, mean(SEAT_SMB_near, 2)+std(SEAT_SMB_near,[],2),'r--')
    plot(siteI_yr, mean(SEAT_SMB_near, 2)-std(SEAT_SMB_near,[],2),'r--')
    
    OIB_SMB_near = OIB_SMB_MC{OIBi_near}(OIB_start(1):...
        OIB_start(1)+length(siteI_yr)-1,:);
    for n = 1:size(OIB_SMB_near, 2)
        h0 = plot(siteI_yr, OIB_SMB_near(:,n), 'm', 'LineWidth', 0.5);
        h0.Color(4) = 0.02;
    end
    h2 = plot(siteI_yr, mean(OIB_SMB_near, 2), 'm', 'LineWidth', 2);
    plot(siteI_yr, mean(OIB_SMB_near, 2)+std(OIB_SMB_near,[],2),'m--')
    plot(siteI_yr, mean(OIB_SMB_near, 2)-std(OIB_SMB_near,[],2),'m--')
    
    core_start = find(core_i.SMB_yr==yr_start);
    core_end = find(core_i.SMB_yr==yr_end);
    coreI_SMB = core_i.SMB(core_start:core_end,:);
    % core4_SMB = movmean(core4_SMB, 3);
    for n = 1:size(coreI_SMB, 2)
        h0 = plot(siteI_yr, coreI_SMB(:,n), 'b', 'LineWidth', 0.5);
        h0.Color(4) = 0.02;
    end
    h3 = plot(siteI_yr, mean(coreI_SMB, 2), 'b', 'LineWidth', 2);
    plot(siteI_yr, mean(coreI_SMB, 2) + std(coreI_SMB, [], 2), 'b--')
    plot(siteI_yr, mean(coreI_SMB, 2) - std(coreI_SMB, [], 2), 'b--')
    xlim([yr_end yr_start])
    xlabel('Calendar Year')
    ylabel('Annual SMB (mm w.e./a)')
    legend([h1 h2 h3], 'SEAT radar', 'OIB radar', 'Firn core')
    hold off
    
%     % Save nearest trace figure
%     f1_name = strcat(cores_loop{i}, '_nearSMB');
%     export_fig(f1, fullfile(out_dir, f1_name), '-png');


    
    f2 =  figure;
    hold on
%     for n = 1:size(SEATi_SMB, 2)
%         h0 = plot(siteI_yr, SEATi_SMB(:,n), 'r', 'LineWidth', 0.5);
%         h0.Color(4) = 0.02;
%     end
%     for n = 1:size(OIBi_SMB, 2)
%         h1 = plot(siteI_yr, OIBi_SMB(:,n), 'm', 'LineWidth', 0.5);
%         h1.Color(4) = 0.02;
%     end
    h1 = plot(siteI_yr, mean(SEATi_SMB, 2), 'r', 'LineWidth', 2);
    plot(siteI_yr, mean(SEATi_SMB, 2) + std(SEATi_SMB, [], 2), 'r--')
    plot(siteI_yr, mean(SEATi_SMB, 2) - std(SEATi_SMB, [], 2), 'r--')
    
    h2 = plot(siteI_yr, mean(OIBi_SMB, 2), 'm', 'LineWidth', 2);
    plot(siteI_yr, mean(OIBi_SMB, 2) + std(OIBi_SMB, [], 2), 'm--')
    plot(siteI_yr, mean(OIBi_SMB, 2) - std(OIBi_SMB, [], 2), 'm--')

    core_start = find(core_i.SMB_yr==yr_start);
    core_end = find(core_i.SMB_yr==yr_end);
    coreI_SMB = core_i.SMB(core_start:core_end,:);
    % core4_SMB = movmean(core4_SMB, 3);
    h3 = plot(siteI_yr, mean(coreI_SMB, 2), 'b', 'LineWidth', 2);
    plot(siteI_yr, mean(coreI_SMB, 2) + std(coreI_SMB, [], 2), 'b--')
    plot(siteI_yr, mean(coreI_SMB, 2) - std(coreI_SMB, [], 2), 'b--')
    xlim([yr_end yr_start])
    xlabel('Calendar Year')
    ylabel('Annual SMB (mm w.e./a)')
    legend([h1 h2 h3], 'SEAT radar', 'OIB radar', 'Firn core')
    hold off
    
%     % Save nearest trace figure
%     f2_name = strcat(cores_loop{i}, '_allSMB');
%     export_fig(f2, fullfile(out_dir, f2_name), '-png');
    
    
    
    


    % Bias relative to the core (nearest trace distribution)
    bias_SEATnear = SEAT_SMB_near - coreI_SMB;
    bias_meanS1 = mean(bias_SEATnear(:));
    bias_stdS1 = std(bias_SEATnear(:));
    bias_semS1 = bias_stdS1/sqrt(numel(bias_SEATnear));
    T_SEATnear = tinv(0.975, numel(bias_SEATnear)-1);
    MoE_SEATnear = T_SEATnear*bias_semS1;
    
    bias_SEATnear = (SEAT_SMB_near - coreI_SMB)./mean(coreI_SMB);
    bias_meanS1_perc = mean(bias_SEATnear(:));
    bias_stdS1_perc = std(bias_SEATnear(:));
    bias_semS1 = bias_stdS1_perc/sqrt(numel(bias_SEATnear));
    T_SEATnear = tinv(0.975, numel(bias_SEATnear)-1);
    MoE_SEATnear_perc = T_SEATnear*bias_semS1;
    
    bias_OIBnear = OIB_SMB_near - coreI_SMB;
    bias_meanO1 = mean(bias_OIBnear(:));
    bias_stdO1 = std(bias_OIBnear(:));
    bias_semO1 = bias_stdO1/sqrt(numel(bias_OIBnear));
    T_OIBnear = tinv(0.975, numel(bias_OIBnear)-1);
    MoE_OIBnear = T_OIBnear*bias_semO1;
    
    bias_OIBnear = (OIB_SMB_near - coreI_SMB)./mean(coreI_SMB);
    bias_meanO1_perc = mean(bias_OIBnear(:));
    bias_stdO1_perc = std(bias_OIBnear(:));
    bias_semO1 = bias_stdO1_perc/sqrt(numel(bias_OIBnear));
    T_OIBnear = tinv(0.975, numel(bias_OIBnear)-1);
    MoE_OIBnear_perc = T_OIBnear*bias_semO1;
    
    % Bias relative to the core (distribution of nearby mean traces)
    bias_SEATi = SEATi_SMB - mean(coreI_SMB,2);
    biasSEATi_mean = mean(bias_SEATi(:));
    biasSEATi_std = std(bias_SEATi(:));
    biasSEATi_sem = biasSEATi_std/sqrt(numel(bias_SEATi));
    T_SEATi = tinv(0.975, numel(bias_SEATi)-1);
    MoE_SEATi = T_SEATi*biasSEATi_sem;
    
    bias_SEATi = (SEATi_SMB - mean(coreI_SMB,2))./mean(coreI_SMB(:));
    biasSEATi_mu_perc = mean(bias_SEATi(:));
    biasSEATi_std_perc = std(bias_SEATi(:));
    biasSEATi_sem = biasSEATi_std_perc/sqrt(numel(bias_SEATi));
    T_SEATi = tinv(0.975, numel(bias_SEATi)-1);
    MoE_SEATi_perc = T_SEATi*biasSEATi_sem;
    
    bias_OIBi = OIBi_SMB - mean(coreI_SMB,2);
    biasOIBi_mean = mean(bias_OIBi(:));
    biasOIBi_std = std(bias_OIBi(:));
    biasOIBi_sem = biasOIBi_std/sqrt(numel(bias_OIBi));
    T_OIBi = tinv(0.975, numel(bias_OIBi)-1);
    MoE_OIBi = T_OIBi*biasOIBi_sem;
    
    bias_OIBi = (OIBi_SMB - mean(coreI_SMB,2))./mean(coreI_SMB(:));
    biasOIBi_mean_perc = mean(bias_OIBi(:));
    biasOIBi_std_perc = std(bias_OIBi(:));
    biasOIBi_sem = biasOIBi_std_perc/sqrt(numel(bias_OIBi));
    T_OIBi = tinv(0.975, numel(bias_OIBi)-1);
    MoE_OIBi_perc = T_OIBi*biasOIBi_sem;
    
    
    bias_mu_abs = [bias_meanS1; biasSEATi_mean; bias_meanO1; biasOIBi_mean];
    bias_mu_perc = [bias_meanS1_perc; biasSEATi_mu_perc; ...
        bias_meanO1_perc; biasOIBi_mean_perc];
    bias_MoE = [MoE_SEATnear; MoE_SEATi; MoE_OIBnear; MoE_OIBi];
    bias_MoE_perc = [MoE_SEATnear_perc; MoE_SEATi_perc; ...
        MoE_OIBnear_perc; MoE_OIBi_perc];
    bias_std = [bias_stdS1; biasSEATi_std; bias_stdO1; biasOIBi_std];
    bias_std_perc = [bias_stdS1_perc; biasSEATi_std_perc; ...
        bias_stdO1_perc; biasOIBi_std_perc];
    stats_i = table(bias_mu_abs, bias_MoE, bias_std, bias_mu_perc,...
        bias_MoE_perc, bias_std_perc, 'VariableNames', ...
        {'mu_abs', 'MoE_abs', 'StdDev_abs',...
        'mu_percent', 'MoE_percent',  'StdDev_percent'},...
        'RowNames', {'SEAT-core (nearest trace)', 'SEAT-core (mean traces)', ...
        'OIB-core (nearest trace)', 'OIB-core (mean traces)'});
    
    
    
    
    % [p] = wanova(coreI_SMB(25,:), SEATi_SMB(25,:), OIBi_SMB(25,:));
    
    %%% CALCULATION OF TRENDS
    
    SEATi_b1 = zeros(1, size(SEATi_SMB, 2));
    SEATi_p = zeros(1, size(SEATi_SMB, 2));
    for j = 1:size(SEATi_SMB, 2)
%         [trend, stats] = robustfit(siteI_yr, SEATi_SMB(:,j));
%         SEATi_b1(j) = trend(2);
%         SEATi_p(j) = stats.p(2);
        [trend,~,~,~,stats] = regress(SEATi_SMB(:,j), ...
            [ones(length(siteI_yr),1) siteI_yr]);
        SEATi_b1(j) = trend(2);
        SEATi_p(j) = stats(3);
    end
    
    OIBi_b1 = zeros(1, size(OIBi_SMB, 2));
    OIBi_p = zeros(1, size(OIBi_SMB, 2));
    for j = 1:size(OIBi_SMB, 2)
%         [trend, stats] = robustfit(siteI_yr, OIBi_SMB(:,j));
%         OIBi_b1(j) = trend(2);
%         OIBi_p(j) = stats.p(2);
        [trend,~,~,~,stats] = regress(OIBi_SMB(:,j), ...
            [ones(length(siteI_yr),1) siteI_yr]);
        OIBi_b1(j) = trend(2);
        OIBi_p(j) = stats(3);
    end
    
    core_trend = zeros(1, size(coreI_SMB,2));
    core_p = zeros(1, size(coreI_SMB,2));
    for j = 1:size(coreI_SMB,2)
%         [trend, stats] = robustfit(siteI_yr, coreI_SMB(:,j));
%         core_trend(j) = trend(2);
%         core_p(j) = stats.p(2);
        [trend_tmp,~,~,~,stats] = regress(coreI_SMB(:,j), ...
            [ones(length(siteI_yr),1) siteI_yr]);
        core_trend(j) = trend_tmp(2);
        core_p(j) = stats(3);
    end
    
    SMB_sites.(cores_loop{i}) = struct('SEAT_SMB', SEATi_SMB, 'SEAT_trend',...
        SEATi_b1, 'SEAT_p', SEATi_p, 'OIB_SMB', OIBi_SMB, 'OIB_trend', ...
        OIBi_b1, 'OIB_p', OIBi_p, 'core_SMB', coreI_SMB, 'core_trend', ...
        core_trend, 'core_p', core_p, 'SMB_yr', siteI_yr, 'bias_stats', ...
        stats_i);
    
%     SEAT_sig = numel(SEATi_p(SEATi_p<0.05))/numel(SEATi_p)
%     OIB_sig = numel(OIBi_p(OIBi_p<0.05))/numel(OIBi_p)
end

% figure
% hold on
% plot(median(SMB_sites.SEAT10_4.SEAT_SMB), SMB_sites.SEAT10_4.SEAT_trend, 'ro')
% plot(median(SMB_sites.SEAT10_4.OIB_SMB), SMB_sites.SEAT10_4.OIB_trend, 'mo')
% plot(median(median(SMB_sites.SEAT10_4.core_SMB)), ...
%     SMB_sites.SEAT10_4.core_trend, 'bo', 'MarkerSize', 10)
% plot(median(SMB_sites.SEAT10_5.SEAT_SMB), SMB_sites.SEAT10_5.SEAT_trend, 'rx')
% plot(median(SMB_sites.SEAT10_5.OIB_SMB), SMB_sites.SEAT10_5.OIB_trend, 'mx')
% plot(median(median(SMB_sites.SEAT10_5.core_SMB)), ...
%     SMB_sites.SEAT10_5.core_trend, 'bx', 'MarkerSize', 10)
% plot(median(SMB_sites.SEAT10_6.SEAT_SMB), SMB_sites.SEAT10_6.SEAT_trend, 'r*')
% plot(median(SMB_sites.SEAT10_6.OIB_SMB), SMB_sites.SEAT10_6.OIB_trend, 'm*')
% plot(median(median(SMB_sites.SEAT10_6.core_SMB)), ...
%     SMB_sites.SEAT10_6.core_trend, 'b*', 'MarkerSize', 10)
% xlabel('Mean annual SMB (mm w.e/a)')
% ylabel('Mean trend in annual SMB (mm/a/a)')
% hold off

f3 = figure('Position', [200 200 1500 800]);
hold on
plot(mean(SMB_sites.SEAT10_4.SEAT_SMB), 100*...
    SMB_sites.SEAT10_4.SEAT_trend./mean(SMB_sites.SEAT10_4.SEAT_SMB),'ro')
plot(mean(SMB_sites.SEAT10_4.OIB_SMB), 100*...
    SMB_sites.SEAT10_4.OIB_trend./mean(SMB_sites.SEAT10_4.OIB_SMB), 'mo')
plot(mean(SMB_sites.SEAT10_4.core_SMB), 100*SMB_sites.SEAT10_4.core_trend...
    ./mean(SMB_sites.SEAT10_4.core_SMB), 'bo')

plot(mean(SMB_sites.SEAT10_5.SEAT_SMB), 100*...
    SMB_sites.SEAT10_5.SEAT_trend./mean(SMB_sites.SEAT10_5.SEAT_SMB),'rx')
plot(mean(SMB_sites.SEAT10_5.OIB_SMB), 100*...
    SMB_sites.SEAT10_5.OIB_trend./mean(SMB_sites.SEAT10_5.OIB_SMB), 'mx')
plot(mean(SMB_sites.SEAT10_5.core_SMB), 100*SMB_sites.SEAT10_5.core_trend...
    ./mean(SMB_sites.SEAT10_5.core_SMB), 'bx')

plot(mean(SMB_sites.SEAT10_6.SEAT_SMB), 100*...
    SMB_sites.SEAT10_6.SEAT_trend./mean(SMB_sites.SEAT10_6.SEAT_SMB),'r*')
plot(mean(SMB_sites.SEAT10_6.OIB_SMB), 100*...
    SMB_sites.SEAT10_6.OIB_trend./mean(SMB_sites.SEAT10_6.OIB_SMB), 'm*')
plot(mean(SMB_sites.SEAT10_6.core_SMB),100*SMB_sites.SEAT10_6.core_trend...
    ./mean(SMB_sites.SEAT10_6.core_SMB), 'b*')
xlabel('Mean annual SMB (mm w.e/a)')
ylabel('Trend in SMB (% change relative to mean SMB)')
hold off

% % Save nearest trace figure
% f3_name = 'trend_comp';
% export_fig(f3, fullfile(out_dir, f3_name), '-png');
% Save SMB-site data for use elsewhere
% save(fullfile(out_dir, 'site_data.mat'), '-struct', 'SMB_sites')

% X = repmat(siteI_yr, 1, size(SEATi_SMB,2));
% x = repmat(siteI_yr, 1, size(SEATi_SMB,2));
% x = reshape(X, numel(X), 1);
% y = reshape(SEATi_SMB, numel(SEATi_SMB), 1);
% [trend,~,~,~,stats] = regress(y, [ones(length(x),1) x]);

%% SEAT and OIB SMB bias tests

% SEATi_N = SEAT_N(SEATi_idx);
% SEATi_E = SEAT_E(SEATi_idx);
% OIBi_N = OIB_N(OIBi_idx);
% OIBi_E = OIB_E(OIBi_idx);

% near_tmp = knnsearch([SEAT_E' SEAT_N'], [OIB_E' OIB_N']);
near_tmp = knnsearch([OIB_E' OIB_N'], [SEAT_E' SEAT_N']);

% [~, idx_O] = unique(near_tmp);
% idx_tmp = 1:length(SEAT_E);
% SEAT_near = idx_tmp(near_tmp(idx_O));
% OIB_near = near_tmp(idx_O);
SEAT_near = 1:length(SEAT_E);
OIB_near = near_tmp;
D_near = diag(pdist2([SEAT_E(SEAT_near)' SEAT_N(SEAT_near)'], ...
    [OIB_E(OIB_near)' OIB_N(OIB_near)']));
D_idx = D_near<=50;
SEAT_near = SEAT_near(D_idx);
OIB_near = OIB_near(D_idx);

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

SEATbias_mean = cell(1,length(bias_yr));
bias_SMB = cell(1,length(bias_yr));
for j=1:length(bias_yr)
    SEATbias_SMB = SEAT_SMB{SEAT_near(j)}(SEAT_start(j):SEAT_end(j));
    SEATbias_mean{j} = mean(SEATbias_SMB);
    OIBbias_SMB = OIB_SMB{OIB_near(j)}(OIB_start(j):OIB_end(j));
    %     bias_SMB{j} = (SEATbias_SMB - OIBbias_SMB)./SEATbias_SMB;
    bias_SMB{j} = (SEATbias_SMB - OIBbias_SMB);
end

% figure
% plot(OIB_E(OIB_near), cellfun(@mean, bias_SMB), '.')


% SMB bias for each individual year in each trace (absolute accumulation)
bias_dist_abs = vertcat(bias_SMB{:});
bias_abs_mu = mean(bias_dist_abs);
bias_abs_std = std(bias_dist_abs);
bias_abs_SEM = bias_abs_std/sqrt(length(bias_dist_abs));
Ts = tinv(0.975, length(bias_dist_abs)-1);
bias_abs_MoE = Ts*bias_abs_SEM;
figure
histogram(bias_dist_abs, 100)

% SMB bias for each individual year in each trace (% change of mean)
bias_dist_perc = cellfun(@(x,y) x/y, bias_SMB, SEATbias_mean, 'UniformOutput', 0);
bias_dist_perc = vertcat(bias_dist_perc{:});
bias_perc_mu = mean(bias_dist_perc);
bias_perc_std = std(bias_dist_perc);
bias_perc_SEM = bias_perc_std/sqrt(length(bias_dist_perc));
Ts = tinv(0.975, length(bias_dist_perc)-1);
bias_perc_MoE = Ts*bias_perc_SEM;
figure
histogram(bias_dist_perc, 100)

bias_mu = [bias_abs_mu; bias_perc_mu];
bias_MoE = [bias_abs_MoE; bias_perc_MoE];
bias_std = [bias_abs_std; bias_perc_std];
bias_stats = table(bias_mu, bias_MoE, bias_std, 'VariableNames', ...
    {'Mean', 'MarginOfError', 'StdDev'}, 'RowNames', ...
    {'SEAT-OIB (absolute SMB)', 'SEAT-OIB (% bias)'});

% save(fullfile(out_dir, 'glob_bias_stats.m'), 'bias_stats')



% % Mean SMB bias for each trace (absolute accumulation)
% bias_dist = cellfun(@mean, bias_SMB);
% bias_bar = mean(bias_dist);
% bias_std = std(bias_dist);
% bias_SEM = bias_std/sqrt(length(bias_dist));
% Ts = tinv([0.025 0.975], length(bias_dist)-1);
% CI = bias_bar + Ts*bias_SEM;
% figure
% histogram(bias_dist, 100)
% 
% % Mean SMB bias for each trace (% change of mean)
% bias_dist = cellfun(@(x,y) mean(x)/y, bias_SMB, SEATbias_mean);
% bias_bar = mean(bias_dist);
% bias_std = std(bias_dist);
% bias_SEM = bias_std/sqrt(length(bias_dist));
% Ts = tinv([0.025 0.975], length(bias_dist)-1);
% CI = bias_bar + Ts*bias_SEM;
% figure
% histogram(bias_dist, 100)


%% Age bias tests

SEAT_ages = [];
for i = 1:length(SEAT_files)
    load(fullfile(SEAT_files(i).folder, SEAT_files(i).name), 'ages');
    SEAT_ages = [SEAT_ages median(ages, 3)];
end


OIB_ages = [];
for i = 1:length(OIB_files)
    load(fullfile(OIB_files(i).folder, OIB_files(i).name), 'ages');
    OIB_ages = [OIB_ages median(ages, 3)];
end
depth = 0:0.02:25;


near_tmp = knnsearch([OIB_E' OIB_N'], [SEAT_E' SEAT_N']);
SEAT_near = 1:length(SEAT_E);
OIB_near = near_tmp;

age_bias = SEAT_ages(:,SEAT_near) - OIB_ages(:,OIB_near);
figure
hold on
plot(depth, mean(age_bias, 2), 'k', 'LineWidth', 2)
plot(depth, mean(age_bias, 2) + std(age_bias, [], 2), 'k--')
plot(depth, mean(age_bias, 2) - std(age_bias, [], 2), 'k--')
hold off

% ERR_tmp = 10*age_bias./(SEAT_ages(1,SEAT_near)-SEAT_ages(:,SEAT_near));
% ERR_decade = mean(ERR_tmp(100:end,:));
res_decade = 10*age_bias(end,:)./(SEAT_ages(1,SEAT_near)-SEAT_ages(end,SEAT_near));

figure
histogram(res_decade, 100)

bias_age_mu = median(res_decade);
bias_age_std = std(res_decade);
bias_SEM = std(bias_age_std)/sqrt(length(res_decade));
Ts = tinv(0.975, length(res_decade)-1);
bias_age_MoE = Ts*bias_SEM;


