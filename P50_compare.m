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
radar_full = load(fullfile(data_path, 'radar/SEAT_Traverses/results_data/OIB_SEAT10_4to10_6_full.mat'));

%% Calculate ages based on different P_50 values (modified from radar_age.m)
tic
% Define horizontal resolution (25 m)
horz_res = 25;

% Calculate continuous layer distances for each layer (accounting for 
% lateral size of stacked radar trace bins)
layers_dist = cellfun(@(x) numel(x)*horz_res, radar_full.layers);

% Map layer prominence-distance values to the location within the radar
% matrix of the ith layer
layer_peaks = zeros(size(radar_full.peaks));
for i = 1:length(radar_full.layers)
    layer_peaks(radar_full.layers{i}) = radar_full.peaks(radar_full.layers{i}).*layers_dist(i);
end

% Define surface age and the year associated with the first pick of the 
% algorithm
age_top = radar_full.collect_date;
yr_pick1 = ceil(radar_full.collect_date - 1);


trace_idx = 1:50:size(layer_peaks, 2);
ages = struct('age500', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age1km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age2km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age3km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age4km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age5km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age6km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
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
    'age18km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age19km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age20km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age21km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age22km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age23km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age24km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
    'age25km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw));
fldnm = fieldnames(ages);
P_50dist = [500 1000:1000:25000];
% ages = struct('age500', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
%     'age1km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
%     'age10km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
%     'age25km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw),...
%     'age50km', zeros(size(layer_peaks, 1), length(trace_idx), Ndraw));
% fldnm = fieldnames(ages);
% P_50dist = [500 1000 10000 25000 50000];

for n = 1:length(fldnm)
    err_out = [];
    for i = 1:length(trace_idx)
        % Assign the 50% likelihood point based on median trace prominence and
        % layer length
        peaks_i = radar_full.peaks(:,trace_idx(i));
        Prom_i = peaks_i(peaks_i>0);
%         P_50 = median(Prom_i)*min([P_50dist(n) 0.5*radar_full.dist(end)]);
        P_50 = min([P_50dist(n) 0.5*radar_full.dist(end)]);
        
        % Assign min/max layer likelihoods, and calculate the logistic rate
        % coefficient
        Po = 0.05;
        K = 1;
        r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
        
        % Get layer prom-distance values and depths for layers in ith trace
        peaks_i = layer_peaks(:,trace_idx(i));
        peaks_idx = peaks_i>0;
        peaks_i = peaks_i(peaks_idx);
        depths_i = radar_full.depth(peaks_idx);
        
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
                ages.(fldnm{n})(:,i,j) = interp1(depths_j, yrs_j, radar_full.depth, 'linear', 'extrap');
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

% ages.age5km = radar_full.ages(:,trace_idx,:);

toc

%% 
clearvars -except radar_full ages* cores Ndraw trace_idx P_50dist
% ages = orderfields(ages, [1 2 3 12 4 5 6 7 8 9 10 11]);
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

end_err = zeros(length(fldnm)-1, length(trace_idx));
for i = 1:length(trace_idx)
    end_err(:,i) = movmean(abs(diff(data_end(:,i))), 2);
%     end_err(:,i) = abs(diff(data_end(:,i)))./abs(diff(P_50))';
end

figure
plot(median(end_err/median(end_err),2))
hold on
plot(median(end_std/median(end_std),2))

% test = (end_err/median(end_err)).*(end_std(2:end,:)/median(end_std(2:end,:)));
% figure
% plot(median(test,2))



% ages = orderfields(ages, [1 2 6 3 4 5]);
% fldnm = fieldnames(ages);
% P_50 = [500 1000 5000 10000 25000 50000];
% C = {'b', 'r', 'm', 'c', 'y', 'g'};
% h = cell2struct(cell(size(fldnm)), fldnm, 1);

% figure
% hold on
% for i = 1:length(fldnm)
%     h.(fldnm{i}) = plot(median(squeeze(ages.(fldnm{i})(end,:,:)), 2), C{i}, ...
%         'LineWidth', 2);
%     plot(median(squeeze(ages.(fldnm{i})(end,:,:)), 2) + ...
%         2*std(squeeze(ages.(fldnm{i})(end,:,:)), [], 2), strcat(C{i}, '--'),...
%         'LineWidth', 0.5)
%     plot(median(squeeze(ages.(fldnm{i})(end,:,:)), 2) - ...
%         2*std(squeeze(ages.(fldnm{i})(end,:,:)), [], 2), strcat(C{i}, '--'),...
%         'LineWidth', 0.5)
% end
% % legend([h.(fldnm{1}) h.(fldnm{2}) h.(fldnm{3}) h.(fldnm{4}) h.(fldnm{5}) h.(fldnm{6})], ...
% %     'P-50=500m', 'P-50=1km', 'P-50=5km', 'P-50=10km', 'P-50=25km', 'P-50=50km')
% hold off