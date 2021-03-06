


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
                %                 data_path = 'C:/Users/durba/Documents/Research/Antarctica/Data/';
                data_path = 'F:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/durba/Documents/MATLAB/Addons/';
        end
        
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile(addon_path, 'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add export_fig to path
addon_struct = dir(fullfile(addon_path, 'altmany-export_fig*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add CReSIS OIB MATLAB reader functions to path
addon_struct = dir(fullfile(addon_path, 'cresis-L1B-matlab-readers*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Define number of Monte Carlo simulations to perform
Ndraw = 100;

horz_res = 25;


%%

% % Name of SEAT core site to generate training data/perform regression
% name = 'SEAT10_6';
% 
% % Load relevant radar data (previously generated using the above section)
% radar = load(fullfile(data_path, 'IceBridge/manual_layers', name, ...
%     strcat('layers_', name, '_opt.mat')));
% 
% % Load manually traced layers for current SEAT core site (generated using
% % 'draw_manual.m')
% tmp = load(fullfile(data_path, 'IceBridge/manual_layers', name, ...
%     strcat('manual_', name, '.mat')));
% man_layers = tmp.man_all;

% Select directory containing PAIPR-results and manual layers (picked using
% the `draw_manual.m` file)
input_dir = uigetdir(fullfile(data_path, "PAIPR-results/vWIP", ...
    "coeffs_branch/manual_layers"), ...
    "Select directory containing input files");

% Load PAIPR radar data
radar = load(fullfile(input_dir, "PAIPR_out"));

% Load manually picked layer indices for above radar echogram
load(fullfile(input_dir, "manual_layers"))
% man_layers = load(fullfile(input_dir, "manual_layers"));
% man_layers = tmp.man_all;

% Define range resolution
range_res = round(mean(diff(radar.depth)),2);

%%

% Find the mean peak prominence (used to scale the prominence-distance
% results)
peak_w = 1/mean(radar.peaks(radar.peaks>0));
dist_w = 1/(size(radar.data_smooth,2)*horz_res);
peak_w = 1;
dist_w = 1;

%Preallocate cell array for auto layer position subscripts
layers_raw = cell(1, max(radar.groups(:)));
for i = 1:length(layers_raw)
    
    % Indices and position subscripts of ith group members
    group_idx = find(radar.groups==i);
    [r,c] = ind2sub(size(radar.data_smooth), group_idx);
    
    % Calculate peak power-distances for each peak in ith group, scaled by
    % the median peak prominence of the radargram
    peak_dist = peak_w*dist_w*radar.peaks(group_idx)*length(c)*horz_res;
%     peak_dist = horz_res*length(c)*radar.peaks(group_idx);
%     d = diff([c(:) r(:)]);
%     peak_dist = horz_res*sum(sqrt(sum(d.*d,2)))*radar.peaks(group_idx);
%     peak_dist = length(c)*radar.peaks(group_idx)/...
%         median(radar.peaks(radar.peaks>0));
    
    % Arrange layer position subscripts (smoothed row values), peak
    % power-distances, and preallocated column for manual layer matches as
    % a matrix, and place in preallocated cell array
    layers_raw{i} = [c movmean(r, 20) peak_dist false(length(c), 1)];
end

% Remove empty cells from layer array
layers_raw = layers_raw(~cellfun(@isempty,layers_raw));

% Sort layers by total power-distance (descending order), so that in
% future calculations operations are performed on brightest/longest 
% layers first
[~, layer_idx] = sort(cellfun(@(x) sum(x(:,3)), layers_raw), 'descend');
layers = layers_raw(layer_idx);

% Allocate array of manual layers with which to search for matching auto
% layers
man_search = man_layers;

% % Remove manual layers that do not extend across at least 50% of the
% % radargram (typically ~15 km)
% max_length = max(cellfun(@length, man_search));
% long_idx = cellfun(@(x) length(x) >= 0.50*max_length, man_search);
% man_search = man_search(long_idx);

for i = 1:length(layers)
    
    % Preallocate matrix for squared sum of errors for distance between
    % manually and auto picked layers
    MAE = zeros(1,length(man_search));
    
    for j = 1:length(man_search)
        
        % Find the col position intersection between ith auto layer and jth
        % manual layer
        [~, ja, jb] = intersect(man_search{j}(:,1), layers{i}(:,1));
        Coord_j = man_search{j}(ja,:);
        layer_j = layers{i}(jb,:);
        
        % Calculate the mean absolute error between ith auto layer and jth
        % manual layer
        MAE(j) = sum(abs(layer_j(:,2) - Coord_j(:,2)))/length(layer_j);
%         SSE(j) = sum((layer_j(:,2)-Coord_j(:,2)).^2)/size(layer_j,1);
    end
    
    % Find the nearest remaining manual layer to the ith auto layer
    [MAE_min, SSE_idx] = min(MAE);
    
    % Determine if nearest manual layer is sufficiently close to ith auto
    % layer
    if MAE_min*range_res <= 0.08
        
        % Remove manual layer nearest to ith auto layer from future
        % match searches
        man_search(SSE_idx) = [];
        
        % Add true label to ith auto layer members column for matched
        % manual layer
        layers{i}(:,4) = true;
    end
    
end

%% Logistic regression

dist_tmp = cellfun(@(x) x(:,3), layers, 'UniformOutput', false);
peak_dist = vertcat(dist_tmp{:});

log_tmp = cellfun(@(x) x(:,4), layers, 'UniformOutput', false);
is_layer = vertcat(log_tmp{:});

[~, dev, stats] = mnrfit(peak_dist, categorical(is_layer));

peaks_plot = (0:0.01:max(peak_dist))';
likelihood = 1./(1 + exp(stats.beta(2)*peaks_plot+stats.beta(1)));

% Diagnostic regression plot
figure
plot(peak_dist, is_layer, '.')
hold on
plot(peaks_plot, likelihood, 'r')


sprintf('%f | %f', stats.beta)


