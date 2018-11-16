%% Why are there such stark differences in trends between adjacent
%  radargrams, even when the mean SMB doesn't have nearly as large as jump?



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

wild = '*.mat';
OIB_files = dir(fullfile(data_path, 'IceBridge/SNO_radar/',...
    '2011/SMB_results/', wild));

radar = load(fullfile(OIB_files(18).folder, OIB_files(18).name));
radar2 = load(fullfile(OIB_files(19).folder, OIB_files(19).name));

yrs = (2010:-1:2010-32+1)';
mSMB = cell2mat(cellfun(@(x) mean(x(1:32,:),2), radar.SMB, ...
    'UniformOutput', false));
mSMB2 = cell2mat(cellfun(@(x) mean(x(1:32,:),2), radar2.SMB, ...
    'UniformOutput', false));

beta = zeros(1,length(mSMB));
rbeta = zeros(1,length(mSMB));
for i=1:length(mSMB)
    coeff = polyfit(yrs, mSMB(:,i), 1);
    beta(i) = coeff(1);
    coeff = robustfit(yrs, mSMB(:,i));
    rbeta(i) = coeff(2);
end

beta2 = zeros(1,length(mSMB2));
rbeta2 = zeros(1,length(mSMB2));
for i=1:length(mSMB2)
    coeff = polyfit(yrs, mSMB2(:,i), 1);
    beta2(i) = coeff(1);
    coeff = robustfit(yrs, mSMB2(:,i));
    rbeta2(i) = coeff(2);
end

E = [radar.Easting radar2.Easting];
data = [radar.data_smooth radar2.data_smooth];



layers = cell(1, max(radar.groups(:)));
for i=1:length(layers)
    [r,c] = find(radar.groups==i);
    layers{i} = [r c];
end
layers = layers(~cellfun(@isempty, layers));

layers2 = cell(1, max(radar2.groups(:)));
for i=1:length(layers2)
    [r,c] = find(radar2.groups==i);
    layers2{i} = [r c];
end
layers2 = layers2(~cellfun(@isempty, layers2));

figure
imagesc(E, radar.depth, data, [-2 2])
hold on
cellfun(@(x) plot(radar.Easting(x(:,2)), radar.depth(x(:,1)), ...
    'r', 'LineWidth', 1.5), layers);
cellfun(@(x) plot(radar2.Easting(x(:,2)), radar2.depth(x(:,1)), ...
    'm', 'LineWidth', 1.5), layers2);

figure
yyaxis left
hold on
plot(radar.Easting, mean(mSMB), 'b')
plot(radar2.Easting, mean(mSMB2), 'c')
yyaxis right
plot(radar.Easting, rbeta, 'r')
plot(radar2.Easting, rbeta2, 'm')
xlim([min(E) max(E)])
hold off

figure
hold on
num_plot = 100;
for n = 1:num_plot
    h0 = plot(yrs, mSMB(:,length(mSMB)-n), 'r', 'LineWidth', 0.5);
    h0.Color(4) = 0.05;
    h0 = plot(yrs, mSMB2(:,length(mSMB2)-n), 'm', 'LineWidth', 0.5);
    h0.Color(4) = 0.05;
end
hold off


%%

% Add gif addon to path
addon_folder = fullfile(addon_path, 'gif_v1.0/');
addpath(genpath(addon_folder))

fig = figure('Position', [10 500 1400 400]);
hold on
scatter(radar.Easting, radar.Northing, 10, mSMB(end,:), 'filled')
scatter(radar2.Easting, radar2.Northing, 10, mSMB2(end,:), 'filled')
scatter(radar.Easting(1), radar.Northing(1), 75, 'kx')
scatter(radar.Easting(end), radar.Northing(end), 75, 'rx')
scatter(radar2.Easting(1), radar2.Northing(1), 75, 'mx')
xlim([min([min(radar.Easting) min(radar2.Easting)]) ...
    max([max(radar.Easting) max(radar2.Easting)])]);
ylim([min([min(radar.Northing) min(radar2.Northing)]) ...
    max([max(radar.Northing) max(radar2.Northing)])]);
colorbar
caxis([150 500])
title(sprintf('%i', yrs(end)))
hold off


gif('test.gif', 'DelayTime', 1, 'LoopCount', 5,'frame', gcf)

for i = length(yrs)-1:-1:1
    fig.NextPlot = 'replace';
    scatter(radar.Easting, radar.Northing, 10, mSMB(i,:), 'filled')
    hold on
    scatter(radar2.Easting, radar2.Northing, 10, mSMB2(i,:), 'filled')
    scatter(radar.Easting(1), radar.Northing(1), 75, 'kx')
    scatter(radar.Easting(end), radar.Northing(end), 75, 'rx')
    scatter(radar2.Easting(1), radar2.Northing(1), 75, 'mx')
    xlim([min([min(radar.Easting) min(radar2.Easting)]) ...
        max([max(radar.Easting) max(radar2.Easting)])]);
    ylim([min([min(radar.Northing) min(radar2.Northing)]) ...
        max([max(radar.Northing) max(radar2.Northing)])]);
    colorbar
    caxis([150 500])
    title(sprintf('%i', yrs(i)))
    hold off
    
    gif
end

