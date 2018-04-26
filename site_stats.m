% This is a script for investigating statistical relationships for radar
% from individual firn core sites

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/Data/';
        addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Addons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Add export_fig to path
addon_folder = strcat(addon_path, 'altmany-export_fig-cafc7c5/');
addpath(genpath(addon_folder))

% Define number of MC realizations
Ndraw = 100;

% Import firn core data
[cores] = import_cores(strcat(data_path, 'Ice-cores/SEAT_cores/', ...
    'DGK_core_data.xlsx'), Ndraw);

out_dir = '/media/durbank/JUMP/Research/Reports-noVC/Figures/';

%%

inputs = {'SEAT10_1' 'SEAT10_3' 'SEAT10_4', 'SEAT10_5', 'SEAT10_6'};

for i = 1:length(inputs)
    file = strcat(data_path, 'radar/SEAT_Traverses/core-site_tests/', ...
        inputs{i}, '_radar.mat');
    radar.(inputs{i}) = load(file);
end

%%

for i = 1:length(inputs)
    
    % Name of ith core site to perform analysis/stats on
    name = inputs{i};
    
    % End year for the shortest record in site dataset
    yr_end = max([min(cores.(name).SMB_yr) cellfun(@min, radar.(name).SMB_yr)]);
    
    % Clip core SMB data to shortest record
    core_end = find(cores.(name).SMB_yr==yr_end);
    cores.(name).SMB_yr = cores.(name).SMB_yr(1:core_end);
    cores.(name).SMB = cores.(name).SMB(1:core_end,:);
    
    % Clip radar SMB data to shortest record
    radar_end = find(radar.(name).SMB_yr{1}==yr_end);
    radar.(name).SMB_yr = cellfun(@(x) x(1:radar_end), radar.(name).SMB_yr,...
        'UniformOutput', 0);
    radar.(name).SMB = cellfun(@(x) x(1:radar_end,:), radar.(name).SMB, ...
        'UniformOutput', 0);
    
    
    
    text_name = strrep(name, '_', '-');
    irand = randi(length(radar.(name).SMB));
    
    f1 = figure('Position', [200 200 1000 550]);
    imagesc(radar.(name).dist, radar.(name).depth, radar.(name).data_smooth, [-2 2])
    colorbar
    xlabel('Distance along profile (m)')
    ylabel('Depth (m)')
    hold on
    plot([radar.(name).dist(irand) radar.(name).dist(irand)], ...
        [0 radar.(name).depth(end)], 'r', 'LineWidth', 2)
    xlim([0 radar.(name).dist(end)])
    ylim([0 radar.(name).depth(end)])
    set(gca, 'Ydir', 'reverse', 'FontSize', 14)
    title(strcat(text_name, ' radargram'));
    hold off
    
    f1_name = strcat(name, '_radargram');
    export_fig(f1, strcat(out_dir, f1_name), '-png');
    close(f1)
    
    f2 = figure('Position', [200 200 750 550]);
    hold on
    for n = 1:length(radar.(name).SMB)
        h0 = plot(radar.(name).depth, median(radar.(name).ages(:,n,:), 3), 'r', 'LineWidth', 0.5);
        h0.Color(4) = 0.01;
    end
    h1 = plot(radar.(name).depth, median(radar.(name).ages(:,irand,:), 3),...
        'r', 'LineWidth', 2);
    plot(radar.(name).depth, median(radar.(name).ages(:,irand,:), 3)...
        + 2*std(squeeze(radar.(name).ages(:,irand,:)), [], 2), 'r--')
    plot(radar.(name).depth, median(radar.(name).ages(:,irand,:), 3)...
        - 2*std(squeeze(radar.(name).ages(:,irand,:)), [], 2), 'r--')
    h2 = plot(cores.(name).depth, mean(cores.(name).ages, 2), 'b');
    plot(cores.(name).depth, mean(cores.(name).ages, 2) + ...
        2*std(cores.(name).ages, [], 2), 'b--')
    plot(cores.(name).depth, mean(cores.(name).ages, 2) - ...
        2*std(cores.(name).ages, [], 2), 'b--')
    title(strcat(text_name, ' age-depth scale'))
    legend([h1 h2],'Sample radar trace', 'Core')
    ylabel('Calendar years')
    xlabel('Depth (m)')
    hold off
    
    f2_name = strcat(name, '_age');
    export_fig(f1, strcat(out_dir, f2_name), '-png');
    close(f2)
    
    
    f3 = figure('Position', [200 200 1300 700]);
    hold on
    title(strcat(text_name, ' annual SMB'))
    for n = 1:length(radar.(name).SMB)
        h0 = plot(radar.(name).SMB_yr{n}, median(radar.(name).SMB{n}, 2), 'r', 'LineWidth', 0.5);
        h0.Color(4) = 0.01;
    end
    h1 = plot(radar.(name).SMB_yr{irand}, median(radar.(name).SMB{irand}, 2),...
        'r', 'LineWidth', 2);
    plot(radar.(name).SMB_yr{irand}, median(radar.(name).SMB{irand}, 2) + ...
        2*std(radar.(name).SMB{irand}, [], 2), 'r--');
    plot(radar.(name).SMB_yr{irand}, median(radar.(name).SMB{irand}, 2) - ...
        2*std(radar.(name).SMB{irand}, [], 2), 'r--');
    h2 = plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2), 'b', 'LineWidth', 2);
    plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) + ...
        2*std(cores.(name).SMB, [], 2), 'b--');
    plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) - ...
        2*std(cores.(name).SMB, [], 2), 'b--');
    legend([h1 h2], 'Sample radar trace', 'Core')
    xlabel('Calendar years')
    ylabel('Annual SMB (mm w.e.)')
    
    f3_name = strcat(name, '_SMB');
    export_fig(f1, strcat(out_dir, f3_name), '-png');
    close(f3)
    
    %%
    % Calculate mean accumulation rate and std at each location
    SMB_mean = cellfun(@mean, radar.(name).SMB, 'UniformOutput', 0);
    SMB_std = cellfun(@std, radar.(name).SMB, 'UniformOutput', 0);
    
    trend = cell(1, length(radar.(name).SMB));
    p_val = cell(1, length(radar.(name).SMB));
    for j = 1:length(radar.(name).SMB)
        
        trend_i = zeros(1, Ndraw);
        p_val_i = zeros(1, Ndraw);
        for k = 1:Ndraw
            % Regression using regress function (average of MC simulations)
            [b, ~, ~, ~, stats] = regress(radar.(name).SMB{j}(:,k), ...
                [ones(length(radar.(name).SMB_yr{j}), 1) radar.(name).SMB_yr{j}]);
            trend_i(k) = b(2);
            p_val_i(k) = stats(3);
        end
        trend{j} = trend_i;
        p_val{j} = p_val_i;
    end
    
    trend_mean = cellfun(@mean, trend);
    trend_std = cellfun(@std, trend);
    
    % Find indices of significant trends (at 95% confidence level)
    p_star = cellfun(@(x) find(x<=0.05), p_val, 'UniformOutput', 0);
    
    % Calculate percentage of MC realizations with significant trends
    p_ratio = cellfun(@(sig, all) length(sig)/length(all), p_star, p_val);
    
    % Regression of SMB trend on mean SMB (all traces)
    [trend_b, trend_int, ~, ~, trend_stats] = regress(trend_mean', ...
        [ones(length(trend_mean), 1) cellfun(@mean, SMB_mean)']);
    [p_b, S_b] = polyfit(cellfun(@mean, SMB_mean), trend_mean, 1);
    [trend_Y, trend_D] = polyconf(p_b, cellfun(@mean, SMB_mean), S_b);
    
    
    
    % Linear trend in SMB for ith core
    trend_core = zeros(1, Ndraw);
    p_core = zeros(1, Ndraw);
    for k = 1:Ndraw
        % Regression using regress function (average of MC simulations)
        [b_core, ~, ~, ~, stats_core] = regress(cores.(name).SMB(:,k), ...
            [ones(length(cores.(name).SMB_yr), 1) cores.(name).SMB_yr]);
        trend_core(k) = b_core(2);
        p_core(k) = stats_core(3);
    end
    
    %% Diagnostic plots for trends
    
    % Mean p-value for each trace
    p_mean = cellfun(@mean, p_val);
    
    % Mean SMB across entire radargram (mean of all realizations for each trace)
    f4 = figure('Position', [200 200 1000 500]);
    hold on
    title(strcat(text_name, ' mean SMB'))
    scatter(radar.(name).Easting, radar.(name).Northing, 50, ...
        cellfun(@mean, SMB_mean), 'filled', 's')
    scatter(cores.(name).Easting, cores.(name).Northing, 150, ...
        mean(mean(cores.(name).SMB, 2)), 'filled')
    hcb = colorbar;
    ylabel(hcb, 'Mean annual SMB (mm/a')
    colormap(cool)
    
    f4_name = strcat(name, '_SMB_map');
    export_fig(f1, strcat(out_dir, f4_name), '-png');
    close(f4)
    
    % Mean trend in SMB across entire radargram
    f5 = figure('Position', [200 200 1000 500]);
    hold on
    title(strcat(text_name, ' mean trend in SMB'))
    scatter(radar.(name).Easting, radar.(name).Northing, 50, trend_mean, ...
        'filled', 's')
    scatter(cores.(name).Easting, cores.(name).Northing, 150, ...
        mean(trend_core), 'filled')
    scatter(radar.(name).Easting(p_mean<0.05), ...
        radar.(name).Northing(p_mean<0.05), 10, 'k+')
    scatter(cores.(name).Easting(mean(p_core)<0.05), ...
        cores.(name).Northing(mean(p_core)<0.05), 10, 'k+')
    hcb = colorbar;
    ylabel(hcb, 'Trend in SMB (mm/a)')
    colormap(cool)
    hold off
    
    f5_name = strcat(name, '_trend_map');
    export_fig(f1, strcat(out_dir, f5_name), '-png');
    close(f5)
    
    
    subset = 1:5:length(SMB_mean);
    
    % Boxplots of SMB mean
    mean_box = (cell2mat(SMB_mean'))';
    
    f6 = figure('Position', [50 550 1800 450]);
    boxplot([(mean(cores.(name).SMB))' mean_box(:,subset)])
    ylabel('Mean SMB (mm/a)')
    set(gca,'Xticklabel',[])
    title(strcat(text_name, ' mean SMB'))
    
    f6_name = strcat(name, '_SMB_box');
    export_fig(f1, strcat(out_dir, f6_name), '-png');
    close(f6)
    
    % Boxplots of trends in SMB
    trend_box = (cell2mat(trend'))';
    
    f7 = figure('Position', [50 50 1800 450]);
    title(strcat(text_name, ' SMB trends'))
    boxplot([trend_core' trend_box(:,subset)])
    ylabel('SMB trend (mm/a)')
    set(gca,'Xticklabel',[])
    title(strcat(text_name, ' SMB trends'))
    
    f7_name = strcat(name, '_trend_box');
    export_fig(f1, strcat(out_dir, f7_name), '-png');
    close(f7)
    
end

%% Core comparisons

% Compare age-depth scales for different cores (currently limited to
% SEAT2010 cores)
f1 = figure('Position', [200 200 750 550]);
hold on
for i = 1:length(cores.name(1:5))
    name_i = cores.name{i};
    plot(cores.(name_i).depth, median(cores.(name_i).ages, 2), 'LineWidth', 2)
end
legend(cores.name(1:5))
xlabel('Depth (m)')
ylabel('Calendar year')
hold off

export_fig(f1, strcat(out_dir, 'core_ages'), '-png');
close(f1)


f2 = figure('Position', [200 200 1300 700]);
hold on
for i = 1:length(cores.name(1:5))
    name_i = cores.name{i};
    plot(cores.(name_i).SMB_yr, median(cores.(name_i).SMB, 2), 'LineWidth', 2)
end
legend(cores.name(1:5))
ylabel('SMB (mm/a)')
xlabel('Calendar year')
hold off

export_fig(f2, strcat(out_dir, 'core_SMB'), '-png');
close(f2)

% f3 = figure('Position', [50 50 1800 450]);
% boxplot([trend_core' trend_box(:,subset)])
% ylabel('SMB trend (mm/a)')
% set(gca,'Xticklabel',[])
% title(strcat(text_name, ' SMB trends'))
%
% export_fig(f3, strcat(out_dir, 'core_trends'), '-png');
% close(f3)