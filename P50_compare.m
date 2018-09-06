% Script to compare the effect of differences in the assigned value of P_50

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

% Load previously processed full OIB radar transect
% radar = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/OIB_SEAT10_4to10_6_full.mat'));
radar = load(fullfile(data_path, ...
    'radar/SEAT_Traverses/results_data/comparisons/dataset_length/OIB_100km.mat'));

%% Calculate ages based on different P_50 values (modified from radar_age.m)

tic
% Define horizontal resolution (25 m)
horz_res = 25;

group_nums = unique(radar.groups(radar.groups>0));
layers = cell(1, length(group_nums));
for i = 1:length(group_nums)
    layers{i} = find(radar.groups == group_nums(i));
end

% Calculate continuous layer distances for each layer (accounting for 
% lateral size of stacked radar trace bins)
layers_dist = cellfun(@(x) numel(x)*horz_res, layers);

% Map layer prominence-distance values to the location within the radar
% matrix of the ith layer
layer_peaks = zeros(size(radar.peaks));
for i = 1:length(layers)
    layer_peaks(layers{i}) = radar.peaks(layers{i}).*layers_dist(i);
end

% Define surface age and the year associated with the first pick of the 
% algorithm
age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);


trace_idx = 1:40:size(layer_peaks, 2);
ages = struct('age1km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age3km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age5km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age7km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age8km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age9km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age10km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age11km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age12km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age13km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age14km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age15km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age16km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age17km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age19km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age21km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age23km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age25km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw));
fldnm = fieldnames(ages);
P_50dist = 1000*[1:2:7 8:17 19:2:25];

for n = 1:length(fldnm)
    err_out = [];
    for i = 1:length(trace_idx)
        % Assign the 50% likelihood point based on median trace prominence and
        % layer length
        peaks_i = radar.peaks(:,trace_idx(i));
        Prom_i = peaks_i(peaks_i>0);
        P_50 = median(Prom_i)*min([P_50dist(n) 0.5*radar.dist(end)]);
        
        % Assign min/max layer likelihoods, and calculate the logistic rate
        % coefficient
        Po = 0.05;
        K = 1;
        r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
        
        % Get layer prom-distance values and depths for layers in ith trace
        peaks_i = layer_peaks(:,trace_idx(i));
        peaks_idx = peaks_i>0;
        peaks_i = peaks_i(peaks_idx);
        depths_i = radar.depth(peaks_idx);
        
        % Likelihood of layer representing a year based on a logistic function
        % with rate (r) calculated above
        likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i));
        
        % Assign MC simulation annual layer presence based on layer likelihood
        % values
        yr_idx = zeros(length(depths_i), Ndraw);
        for j = 1:length(depths_i)
            R = rand(Ndraw, 1) <= likelihood(j);
            yr_idx(j,:) = R;
        end
        
        for j = 1:Ndraw
            depths_j = [0; depths_i(logical(yr_idx(:,j)))];
            yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
            try
                ages.(fldnm{n})(:,i,j) = interp1(depths_j, yrs_j, radar.depth, 'linear', 'extrap');
            catch
                sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
                err_out = [err_out j];
            end
        end
        if ~isempty(err_out)
            ages.(fldnm{n})(:,i,err_out) = repmat(sum(squeeze(ages_500(:,i,:)), 2)./...
                sum(squeeze(ages_500(:,i,:))~=0, 2), 1, length(err_out));
        end
    end
end
toc

%% 
clearvars -except radar ages* cores Ndraw trace_idx P_50dist
fldnm = fieldnames(ages);
P_50 = P_50dist;
h = cell2struct(cell(size(fldnm)), fldnm, 1);

figure
hold on
data_end = zeros(length(fldnm), length(trace_idx));
end_std = zeros(length(fldnm), length(trace_idx));
for i = 1:length(fldnm)
    data_end(i,:) = median(squeeze(ages.(fldnm{i})(end,:,:)), 2);
    end_std(i,:) = std(squeeze(ages.(fldnm{i})(end,:,:)), [], 2);
    h.(fldnm{i}) = plot(median(squeeze(ages.(fldnm{i})(end,:,:)), 2), 'LineWidth', 2);
end
legend(fldnm)
hold off

end_err = zeros(length(fldnm), length(trace_idx));
for i = 1:length(trace_idx)
    end_err(1,i) = abs(diff(data_end(1:2,i)));
    for j = 2:length(fldnm)-1
        end_err(j,i) = mean(abs(diff(data_end(j-1:j+1,i))));
    end
    end_err(end,i) = abs(diff(data_end(end-1:end,i)));
end

figure
hold on
plot(P_50, median(end_err,2), 'k')
plot(P_50, median(end_err(:,1:20),2), 'b')
plot(P_50, median(end_err(:,57:77),2), 'r')
plot(P_50, median(end_err(:,end-20:end),2), 'm')
legend('all', 'SEAT10-4', 'SEAT10-5', 'SEAT10-6')
hold off
% figure
% plot(P_50(2:end), median(end_err,2))
% hold on
% plot(P_50, median(end_std,2))
% test = median(end_err,2);
% plot(P_50, [test(1); test] + median(end_std,2))

%%

% Compare results to cores
sites = {'SEAT10_4' 'SEAT10_5' 'SEAT10_6'};
for i = 1:length(sites)
    name = sites{i};
    dist_radar = pdist2([cores.(name).Easting  cores.(name).Northing], ...
        [radar.Easting(trace_idx)', radar.Northing(trace_idx)']);
    [~, radar_near] = min(dist_radar);
    figure
    hold on
    h1 = plot(cores.(name).depth, mean(cores.(name).ages, 2), 'k', 'LineWidth', 2);
    plot(cores.(name).depth, mean(cores.(name).ages, 2) + ...
        2*std(cores.(name).ages, [], 2), 'k--')
    plot(cores.(name).depth, mean(cores.(name).ages, 2) - ...
        2*std(cores.(name).ages, [], 2), 'k--')
    h2 = plot(radar.depth, squeeze(median(ages.age8km(:,radar_near,:),3)), 'b', 'LineWidth', 2);
    plot(radar.depth, squeeze(median(ages.age8km(:,radar_near,:),3)) + ...
        2*std(ages.age8km(:,radar_near,:), [], 3), 'b--')
    plot(radar.depth, squeeze(median(ages.age8km(:,radar_near,:),3)) - ...
        2*std(ages.age8km(:,radar_near,:), [], 3), 'b--')
    h3 = plot(radar.depth, squeeze(median(ages.age5km(:,radar_near,:),3)), 'c');
    plot(radar.depth, squeeze(median(ages.age5km(:,radar_near,:),3)) + ...
        2*std(ages.age5km(:,radar_near,:), [], 3), 'c--')
    plot(radar.depth, squeeze(median(ages.age5km(:,radar_near,:),3)) - ...
        2*std(ages.age5km(:,radar_near,:), [], 3), 'c--')
    h4 = plot(radar.depth, squeeze(median(ages.age10km(:,radar_near,:),3)), 'r');
    h5 = plot(radar.depth, squeeze(median(ages.age15km(:,radar_near,:),3)), 'm');
    legend([h1 h2 h3 h4 h5], 'core', 'P50=8km', 'P50=5km', 'P50=10km', 'P50=15km')
    xlabel('Depth (m)')
    ylabel('Calendar year')
    hold off
end

