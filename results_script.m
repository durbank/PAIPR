% Script for assessing and comparing results of accum-radar, and generating
% figures for the methods paper

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
OIB_SMB = cellfun(@(x) median(x, 2), OIB_SMB_MC, 'UniformOutput', 0);
OIB_std = cellfun(@(x) std(x, [], 2), OIB_SMB_MC, 'UniformOutput', 0);

% %% Compare individual SMB distributions within 1 km of site 4 (core site SEAT10-4)

% near_dist = 100;
%
% D_SEATi = pdist2([cores.SEAT10_4.Easting, cores.SEAT10_4.Northing],...
%     [SEAT_E' SEAT_N']);
% SEATi_idx = D_SEATi <= near_dist;
%
% D_OIBi = pdist2([cores.SEAT10_4.Easting, cores.SEAT10_4.Northing],...
%     [OIB_E' OIB_N']);
% OIBi_idx = D_OIBi <= near_dist;
%
% SEAT4_dist = SEAT_SMB_MC(SEATi_idx);
% SEAT4_yr = SEAT_yr(SEATi_idx);


%% Compare median SMB estimates within 10 km of site 4 (core site SEAT10-4)

cores_loop = {'SEAT10_4' 'SEAT10_5' 'SEAT10_6'};
near_dist = 10000;
SMB_sites = struct();

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
    
    figure
    hold on
    for n = 1:size(SEATi_SMB, 2)
        h0 = plot(siteI_yr, SEATi_SMB(:,n), 'r', 'LineWidth', 0.5);
        h0.Color(4) = 0.02;
    end
    for n = 1:size(OIBi_SMB, 2)
        h1 = plot(siteI_yr, OIBi_SMB(:,n), 'm', 'LineWidth', 0.5);
        h1.Color(4) = 0.02;
    end
    h2 = plot(siteI_yr, median(SEATi_SMB, 2), 'r', 'LineWidth', 2);
    SEAT_SMB_near = SEAT_SMB{SEATi_near}(SEAT_start(1):SEAT_start(1)+length(siteI_yr)-1);
    h3 = plot(siteI_yr, SEAT_SMB_near, 'r');
    plot(siteI_yr, SEAT_SMB_near + SEAT_std{SEATi_near}...
        (SEAT_start(1):SEAT_start(1)+length(siteI_yr)-1), 'r--', 'LineWidth', 0.5);
    plot(siteI_yr, SEAT_SMB_near - SEAT_std{SEATi_near}...
        (SEAT_start(1):SEAT_start(1)+length(siteI_yr)-1), 'r--', 'LineWidth', 0.5);
    
    h4 = plot(siteI_yr, median(OIBi_SMB, 2), 'm', 'LineWidth', 2);
    OIB_SMB_near = OIB_SMB{OIBi_near}(OIB_start(1):OIB_start(1)+length(siteI_yr)-1);
    h3 = plot(siteI_yr, OIB_SMB_near, 'm');
    plot(siteI_yr, OIB_SMB_near + OIB_std{OIBi_near}...
        (OIB_start(1):OIB_start(1)+length(siteI_yr)-1), 'm--', 'LineWidth', 0.5);
    plot(siteI_yr, OIB_SMB_near - OIB_std{OIBi_near}...
        (OIB_start(1):OIB_start(1)+length(siteI_yr)-1), 'm--', 'LineWidth', 0.5);
    
    core_start = find(core_i.SMB_yr==yr_start);
    core_end = find(core_i.SMB_yr==yr_end);
    coreI_SMB = core_i.SMB(core_start:core_end,:);
    % core4_SMB = movmean(core4_SMB, 3);
    h = plot(siteI_yr, median(coreI_SMB, 2), 'b', 'LineWidth', 1);
    plot(siteI_yr, median(coreI_SMB, 2) + std(coreI_SMB, [], 2), 'b--', 'LineWidth', 0.5);
    plot(siteI_yr, median(coreI_SMB, 2) - std(coreI_SMB, [], 2), 'b--', 'LineWidth', 0.5);
    xlim([yr_end yr_start])
    hold off
    
    % [p] = wanova(coreI_SMB(25,:), SEATi_SMB(25,:), OIBi_SMB(25,:));
    
    %%% CALCULATION OF TRENDS
    
    SEATi_b1 = zeros(1, size(SEATi_SMB, 2));
    SEATi_p = zeros(1, size(SEATi_SMB, 2));
    
    for j = 1:size(SEATi_SMB, 2)
        %     [trend, stats] = robustfit(siteI_yr, SEATi_SMB(:,j));
        %     SEATi_b1(j) = trend(2);
        %     SEATi_p(j) = stats.p(2);
        [trend,~,~,~,stats] = regress(SEATi_SMB(:,j), [ones(length(siteI_yr),1) siteI_yr]);
        SEATi_b1(j) = trend(2);
        SEATi_p(j) = stats(3);
    end
    
    OIBi_b1 = zeros(1, size(OIBi_SMB, 2));
    OIBi_p = zeros(1, size(OIBi_SMB, 2));
    
    for j = 1:size(OIBi_SMB, 2)
        [trend,~,~,~,stats] = regress(OIBi_SMB(:,j), [ones(length(siteI_yr),1) siteI_yr]);
        OIBi_b1(j) = trend(2);
        OIBi_p(j) = stats(3);
    end
    
    
%     [trend_tmp,~,~,~,stats] = regress(median(coreI_SMB,2), [ones(length(siteI_yr),1) siteI_yr]);
%     core_trend = trend_tmp(2);
%     core_p = stats(3)
core_trend = zeros(1, size(coreI_SMB,2));
core_p = zeros(1, size(coreI_SMB,2));
for j = 1:size(coreI_SMB,2)
    [trend_tmp,~,~,~,stats] = regress(coreI_SMB(:,j), ...
        [ones(length(siteI_yr),1) siteI_yr]);
    core_trend(j) = trend_tmp(2);
    core_p(j) = stats(3);
end

    SMB_sites.(cores_loop{i}) = struct('SEAT_SMB', SEATi_SMB, 'SEAT_trend',...
        SEATi_b1, 'SEAT_p', SEATi_p, 'OIB_SMB', OIBi_SMB, 'OIB_trend', ...
        OIBi_b1, 'OIB_p', OIBi_p, 'core_SMB', coreI_SMB, 'core_trend', ...
        core_trend, 'core_p', core_p, 'SMB_yr', siteI_yr);
    
    SEAT_sig = numel(SEATi_p(SEATi_p<0.05))/numel(SEATi_p)
    OIB_sig = numel(OIBi_p(OIBi_p<0.05))/numel(OIBi_p)
    
%     figure
%     hold on
%     plot(median(SEATi_SMB), SEATi_b1./median(SEATi_SMB), 'r.')
%     plot(median(OIBi_SMB), OIBi_b1./median(OIBi_SMB), 'm.')
%     plot(median(median(coreI_SMB,2)), core_trend/median(median(coreI_SMB,2)), 'bx')
%     hold off
    
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

figure
hold on
plot(median(SMB_sites.SEAT10_4.SEAT_SMB), 100*...
    SMB_sites.SEAT10_4.SEAT_trend./median(SMB_sites.SEAT10_4.SEAT_SMB),'ro')
plot(median(SMB_sites.SEAT10_4.OIB_SMB), 100*...
    SMB_sites.SEAT10_4.OIB_trend./median(SMB_sites.SEAT10_4.OIB_SMB), 'mo')
plot(median(SMB_sites.SEAT10_4.core_SMB), 100*SMB_sites.SEAT10_4.core_trend...
    ./median(SMB_sites.SEAT10_4.core_SMB), 'bo')

plot(median(SMB_sites.SEAT10_5.SEAT_SMB), 100*...
    SMB_sites.SEAT10_5.SEAT_trend./median(SMB_sites.SEAT10_5.SEAT_SMB),'rx')
plot(median(SMB_sites.SEAT10_5.OIB_SMB), 100*...
    SMB_sites.SEAT10_5.OIB_trend./median(SMB_sites.SEAT10_5.OIB_SMB), 'mx')
plot(median(SMB_sites.SEAT10_5.core_SMB), 100*SMB_sites.SEAT10_5.core_trend...
    ./median(SMB_sites.SEAT10_5.core_SMB), 'bx')

plot(median(SMB_sites.SEAT10_6.SEAT_SMB), 100*...
    SMB_sites.SEAT10_6.SEAT_trend./median(SMB_sites.SEAT10_6.SEAT_SMB),'r*')
plot(median(SMB_sites.SEAT10_6.OIB_SMB), 100*...
    SMB_sites.SEAT10_6.OIB_trend./median(SMB_sites.SEAT10_6.OIB_SMB), 'm*')
plot(median(SMB_sites.SEAT10_6.core_SMB),100*SMB_sites.SEAT10_6.core_trend...
    ./median(SMB_sites.SEAT10_6.core_SMB), 'b*')
xlabel('Mean annual SMB (mm w.e/a)')
ylabel('Trend in SMB (% change relative to mean SMB)')
hold off

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

SEATbias_med = cell(1,length(bias_yr));
bias_SMB = cell(1,length(bias_yr));
for j=1:length(bias_yr)
    SEATbias_SMB = SEAT_SMB{SEAT_near(j)}(SEAT_start(j):SEAT_end(j));
    SEATbias_med{j} = median(SEATbias_SMB);
    OIBbias_SMB = OIB_SMB{OIB_near(j)}(OIB_start(j):OIB_end(j));
    %     bias_SMB{j} = (SEATbias_SMB - OIBbias_SMB)./SEATbias_SMB;
    bias_SMB{j} = (SEATbias_SMB - OIBbias_SMB);
end

figure
histogram(vertcat(bias_SMB{:}), 100)

figure
plot(OIB_E(OIB_near), cellfun(@mean, bias_SMB), '.')

bias = mean(vertcat(bias_SMB{:}));




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

bias_med = median(res_decade);
bias_SEM = std(res_decade)/sqrt(length(res_decade));
Ts = tinv([0.025 0.795], length(res_decade)-1);
CI = bias_med + Ts*bias_SEM;


