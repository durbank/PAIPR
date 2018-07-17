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

%%

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

Ndraw = 100;

% Load previously processed radar accumulation data
% radar = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/SEAT10_4to10_6_processed.mat'));
radar_OIB = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/OIB_SEAT10_4to10_6.mat'));


% % Load individual segment data, and other pertinent data
% seg1 = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/SEAT10_4to10_6_seg1.mat'));
% seg2 = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/SEAT10_4to10_6_seg2.mat'));
% seg3 = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/SEAT10_4to10_6_seg3.mat'));
% Ndraw = seg1.Ndraw;
% cores = seg1.cores;
% 
% % Combine radar segments into a single radar structure, clipping the ends
% % off of each segment
% clip = 100;
% radar = struct('collect_date', seg1.radar.collect_date, 'Easting', ...
%     [seg1.radar.Easting(clip:end-clip) seg2.radar.Easting(clip:end-clip) seg3.radar.Easting(clip:end-clip)],...
%     'Northing', [seg1.radar.Northing(clip:end-clip) seg2.radar.Northing(clip:end-clip) seg3.radar.Northing(clip:end-clip)],...
%     'dist', [seg1.radar.dist(clip:end-clip) seg2.radar.dist(clip:end-clip) seg3.radar.dist(clip:end-clip)],...
%     'depth', seg1.radar.depth, ...
%     'data_smooth', [seg1.radar.data_smooth(:,clip:end-clip) seg2.radar.data_smooth(:,clip:end-clip) seg3.radar.data_smooth(:,clip:end-clip)],...
%     'likelihood', [seg1.radar.likelihood(:,clip:end-clip) seg2.radar.likelihood(:,clip:end-clip) seg3.radar.likelihood(:,clip:end-clip)],...
%     'ages', [seg1.radar.ages(:,clip:end-clip,:) seg2.radar.ages(:,clip:end-clip,:) seg3.radar.ages(:,clip:end-clip,:)]);
% radar.SMB_yr = [seg1.radar.SMB_yr(clip:end-clip) seg2.radar.SMB_yr(clip:end-clip) seg3.radar.SMB_yr(clip:end-clip)];
% radar.SMB = [seg1.radar.SMB(clip:end-clip) seg2.radar.SMB(clip:end-clip) seg3.radar.SMB(clip:end-clip)];

%% Index map and study site

labels = strrep(cores.name, '_', '-');
basins = shaperead(strcat(data_path, 'DEMs/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'));

Easting_lims = [min(cores.Easting)-0.10*(max(cores.Easting)-min(cores.Easting)) ...
    max(cores.Easting)+0.15*(max(cores.Easting)-min(cores.Easting))];
Northing_lims = [min(cores.Northing)-0.1*(max(cores.Northing)-min(cores.Northing)) ...
    max(cores.Northing)+0.1*(max(cores.Northing)-min(cores.Northing))];
elev = cryosat2_data(Easting_lims, Northing_lims);

figure('Position', [10 10 1400 800])
hold on
h0 = image(Easting_lims, Northing_lims, elev, 'CDataMapping', 'scaled');
colormap(gray)
mapshow(basins, 'FaceAlpha', 0, 'LineWidth', 3)
h1 = scatter(cores.Easting, cores.Northing, 100, 'b', 'filled');
h2 = plot(radar_OIB.Easting(1), radar_OIB.Northing(1), 'r', 'LineWidth', 2);     % Correctly display radar as line in legend
plot(radar_OIB.Easting, radar_OIB.Northing, 'm.', 'MarkerSize', 0.10)
text(cores.Easting, cores.Northing, strcat('\leftarrow', labels), ...
    'FontSize', 18, 'Interpreter', 'tex');
h3 = plot(cores.Easting(1), cores.Northing(1), 'k', 'LineWidth', 2);
% h4 = plot(OIB.Easting(1), OIB.Northing(1), 'm', 'LineWidth', 2);
% plot(OIB.Easting, OIB.Northing, 'm.', 'MarkerSize', 0.05)
c0 = colorbar;
c0.Label.String = 'Elevation (m asl)';
c0.Label.FontSize = 18;
graticuleps(-81:0.5:-77,-125:2:-105, 'c')
xlim(Easting_lims)
ylim(Northing_lims)
scalebarps
box on
mapzoomps('ne', 'insetsize', 0.30)
legend([h1 h2 h3], 'Firn cores', 'Radar transects', 'WAIS Divide', 'Location', 'northwest')
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 18)
hold off


%% Core site comparisons

inputs = {'SEAT10_4', 'SEAT10_5', 'SEAT10_6'};

for i = 1:length(inputs)
    
    % Name of ith core site to perform analysis/stats on
    name = inputs{i};
    file = fullfile(data_path, 'radar/SEAT_Traverses/results_data', ...
        strcat('grid', name, '.mat'));
    radar_i = load(file);
    
    dist_OIB = pdist2([ cores.(name).Easting  cores.(name).Northing], ...
        [radar_OIB.Easting', radar_OIB.Northing']);
    trace_idx = 1:length(radar_OIB.SMB);
    OIB_idx = trace_idx(dist_OIB<=5000);
    [~, OIB_near] = min(dist_OIB);
    %     core_pos = radar.dist(near_idx);
    
    dist_SEAT = pdist2([ cores.(name).Easting  cores.(name).Northing], ...
        [radar_i.Easting', radar_i.Northing']);
%     trace_idx = 1:length(radar_OIB.SMB);
%     SEAT_idx = trace_idx(dist_SEAT<=5000);
    [~, SEAT_near] = min(dist_SEAT);
    
    text_name = strrep(name, '_', '-');
    f1 = figure('Position', [200 200 1000 550]);
    imagesc(radar_OIB.Easting(OIB_idx), radar_OIB.depth, radar_OIB.data_smooth(:,OIB_idx), [-2 2])
    colorbar
    xlabel('Distance along profile (m)')
    ylabel('Depth (m)')
    hold on
    plot([cores.(name).Easting cores.(name).Easting], ...
        [radar_OIB.depth(1) radar_OIB.depth(end)], 'b', 'LineWidth', 2)
    set(gca, 'Ydir', 'reverse', 'FontSize', 14)
    title(strcat(text_name, ' radargram'));
    hold off
    f1_name = strcat(name, '_radargram');
%     export_fig(f1, strcat(out_dir, f1_name), '-png');
%     close(f1)

    f2 = figure('Position', [200 200 750 550]);
    hold on
    for n = OIB_idx
        h0 = plot(radar_OIB.depth, median(radar_OIB.ages(:,n,:), 3), 'm', 'LineWidth', 0.5);
        h0.Color(4) = 0.02;
    end
    for n = 1:length(radar_i.SMB)
        h0 = plot(radar_i.depth, median(radar_i.ages(:,n,:), 3), 'r', 'LineWidth', 0.5);
        h0.Color(4) = 0.02;
    end
    h1 = plot(radar_OIB.depth, median(median(radar_OIB.ages(:,OIB_idx,:), 3), 2),...
        'm--', 'LineWidth', 2);
    h2 = plot(radar_i.depth, median(median(radar_i.ages, 3), 2),...
        'r--', 'LineWidth', 2);
    h3 = plot(cores.(name).depth, mean(cores.(name).ages, 2), 'b');
    plot(cores.(name).depth, mean(cores.(name).ages, 2) + ...
        2*std(cores.(name).ages, [], 2), 'b--')
    plot(cores.(name).depth, mean(cores.(name).ages, 2) - ...
        2*std(cores.(name).ages, [], 2), 'b--')
    title(strcat(text_name, ' age-depth scale'))
    legend([h1 h2 h3],'SEAT traces', 'OIB traces', 'Core')
    ylabel('Calendar years')
    xlabel('Depth (m)')
    hold off
    f2_name = strcat(name, '_age');
%     export_fig(f1, strcat(out_dir, f2_name), '-png');
%     close(f2)
    
    figure('Position', [200 200 750 550]);
    hold on
%     for n = 1:Ndraw
%         h0 = plot(radar_OIB.depth, radar_OIB.ages(:,OIB_near,n), 'm', 'LineWidth', 0.5);
%         h0.Color(4) = 0.02;
%     end
%     for n = 1:Ndraw
%         h0 = plot(radar_i.depth, radar_i.ages(:,SEAT_near,n), 'r', 'LineWidth', 0.5);
%         h0.Color(4) = 0.02;
%     end
    h1 = plot(radar_OIB.depth, median(radar_OIB.ages(:,OIB_near,:), 3),...
        'm', 'LineWidth', 2);
    plot(radar_OIB.depth, median(radar_OIB.ages(:,OIB_near,:), 3)...
        + 2*std(squeeze(radar_OIB.ages(:,OIB_near,:)), [], 2), 'm--')
    plot(radar_OIB.depth, median(radar_OIB.ages(:,OIB_near,:), 3)...
        - 2*std(squeeze(radar_OIB.ages(:,OIB_near,:)), [], 2), 'm--')
    h2 = plot(radar_i.depth, median(radar_i.ages(:,SEAT_near,:), 3),...
        'r', 'LineWidth', 2);
    plot(radar_i.depth, median(radar_i.ages(:,SEAT_near,:), 3)...
        + 2*std(squeeze(radar_i.ages(:,SEAT_near,:)), [], 2), 'r--')
    plot(radar_i.depth, median(radar_i.ages(:,SEAT_near,:), 3)...
        - 2*std(squeeze(radar_i.ages(:,SEAT_near,:)), [], 2), 'r--')
    h3 = plot(cores.(name).depth, mean(cores.(name).ages, 2), 'b');
    plot(cores.(name).depth, mean(cores.(name).ages, 2) + ...
        2*std(cores.(name).ages, [], 2), 'b--')
    plot(cores.(name).depth, mean(cores.(name).ages, 2) - ...
        2*std(cores.(name).ages, [], 2), 'b--')
    title(strcat(text_name, ' age-depth scale'))
    legend([h1 h2 h3],'OIB traces', 'SEAT traces', 'Core')
    ylabel('Calendar years')
    xlabel('Depth (m)')
    hold off

    f3 = figure('Position', [200 200 1300 700]);
    hold on
    title(strcat(text_name, ' annual SMB'))
    for n = OIB_idx
        h0 = plot(radar_OIB.SMB_yr{n}, median(radar_OIB.SMB{n}, 2), 'm', 'LineWidth', 0.5);
        h0.Color(4) = 0.02;
    end
    for n = 1:length(radar_i.SMB)
        h0 = plot(radar_i.SMB_yr{n}, median(radar_i.SMB{n}, 2), 'r', 'LineWidth', 0.5);
        h0.Color(4) = 0.02;
    end
    yr_end = min(cellfun(@length, radar_OIB.SMB_yr(OIB_idx)));
    SMB_data = cell2mat(cellfun(@(x) x(1:yr_end), radar_OIB.SMB(OIB_idx), 'Uniform', false)')';
    h1 = plot(radar_OIB.SMB_yr{OIB_idx(1)}(1:yr_end), median(SMB_data, 2), 'm--', 'LineWidth', 2);
    yr_end = min(cellfun(@length, radar_i.SMB_yr));
    SMB_data = cell2mat(cellfun(@(x) x(1:yr_end), radar_i.SMB, 'Uniform', false)')';
    h2 = plot(radar_i.SMB_yr{1}(1:yr_end), median(SMB_data, 2), 'r--', 'LineWidth', 2);
    h3 = plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2), 'b', 'LineWidth', 2);
    plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) + ...
        2*std(cores.(name).SMB, [], 2), 'b--');
    plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) - ...
        2*std(cores.(name).SMB, [], 2), 'b--');
    legend([h1 h2 h3], 'OIB traces', 'SEAT traces', 'Core')
    xlabel('Calendar years')
    ylabel('Annual SMB (mm w.e.)')
    f3_name = strcat(name, '_SMB');
%     export_fig(f1, strcat(out_dir, f3_name), '-png');
%     close(f3)

    figure('Position', [200 200 1300 700]);
    hold on
    title(strcat(text_name, ' annual SMB'))
%     for n = 1:Ndraw
%         h0 = plot(radar_OIB.SMB_yr{OIB_near}, radar_OIB.SMB{OIB_near}(:,n), 'm', 'LineWidth', 0.5);
%         h0.Color(4) = 0.02;
%     end
%     for n = 1:Ndraw
%         h0 = plot(radar_i.SMB_yr{SEAT_near}, radar_i.SMB{SEAT_near}(:,n), 'r', 'LineWidth', 0.5);
%         h0.Color(4) = 0.02;
%     end
    h1 = plot(radar_OIB.SMB_yr{OIB_near}, median(radar_OIB.SMB{OIB_near}, 2), 'm', 'LineWidth', 2);
    plot(radar_OIB.SMB_yr{OIB_near}, median(radar_OIB.SMB{OIB_near}, 2) + ...
        2*std(radar_OIB.SMB{OIB_near}, [], 2), 'm--');
    plot(radar_OIB.SMB_yr{OIB_near}, median(radar_OIB.SMB{OIB_near}, 2) - ...
        2*std(radar_OIB.SMB{OIB_near}, [], 2), 'm--');
    h2 = plot(radar_i.SMB_yr{SEAT_near}, median(radar_i.SMB{SEAT_near}, 2), 'r', 'LineWidth', 2);
    plot(radar_i.SMB_yr{SEAT_near}, median(radar_i.SMB{SEAT_near}, 2) + ...
        2*std(radar_i.SMB{SEAT_near}, [], 2), 'r--');
    plot(radar_i.SMB_yr{SEAT_near}, median(radar_i.SMB{SEAT_near}, 2) - ...
        2*std(radar_i.SMB{SEAT_near}, [], 2), 'r--');
    h3 = plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2), 'b', 'LineWidth', 2);
    plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) + ...
        2*std(cores.(name).SMB, [], 2), 'b--');
    plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2) - ...
        2*std(cores.(name).SMB, [], 2), 'b--');
    legend([h1 h2 h3], 'OIB traces', 'SEAT traces', 'Core')
    xlabel('Calendar years')
    ylabel('Annual SMB (mm w.e.)')
    hold off

end
%%
radar = radar_OIB;

% Calculate mean accumulation rate and std at each location
SMB_mean = cellfun(@mean, radar.SMB, 'UniformOutput', 0);
SMB_std = cellfun(@std, radar.SMB, 'UniformOutput', 0);

trend = cell(1, length(radar.SMB));
p_val = cell(1, length(radar.SMB));
for i = 1:length(radar.SMB)
    
    trend_i = zeros(1, Ndraw);
    p_val_i = zeros(1, Ndraw);
    for j = 1:Ndraw
        % Regression using regress function (average of MC simulations)
        [b, ~, ~, ~, stats] = regress(radar.SMB{i}(:,j), ...
            [ones(length(radar.SMB_yr{i}), 1) radar.SMB_yr{i}]);
        trend_i(j) = b(2);
        p_val_i(j) = stats(3);
    end
    trend{i} = trend_i;
    p_val{i} = p_val_i;
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

%%

subset = 1:length(SMB_mean);

% Boxplots of SMB mean
mean_box = (cell2mat(SMB_mean'))';

f6 = figure('Position', [50 550 1800 450]);
% boxplot([(mean(cores.(name).SMB))' mean_box(:,subset)])
boxplot(mean_box(:,subset))
ylabel('Mean SMB (mm/a)')
set(gca,'Xticklabel',[])
title('Mean SMB')

f6_name = 'SMB_box';
% export_fig(f1, strcat(out_dir, f6_name), '-png');
% close(f6)

% Boxplots of trends in SMB
trend_box = (cell2mat(trend'))';

f7 = figure('Position', [50 50 1800 450]);
title('SMB trends')
% boxplot([trend_core' trend_box(:,subset)])
boxplot(trend_box(:,subset))
ylabel('SMB trend (mm/a)')
set(gca,'Xticklabel',[])
title('SMB trends')

f7_name = strcat(name, '_trend_box');
% export_fig(f1, strcat(out_dir, f7_name), '-png');
% close(f7)