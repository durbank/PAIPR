
tic

% Create cell arrays of peak prominences and corresponding depths. Then
% loop through the columns (i=2:size(cell,2)-1), where for each peak depth,
% compare the nearest depths on both sides of that column, for any peak
% depths that overlap within a tolerance (1/2*bin_res), set those depths to
% the depth of the strongest prominence.

load('/Users/Durbank/Documents/MATLAB/Research/tmp_data/RT-tests.mat')

%% Find the depth+prominence of peaks within each smoothed radar trace

% Preallocate arrays for various components
peaks = zeros(size(radar.data_smooth));
peak_width = zeros(size(radar.data_smooth));
Proms = cell(1, size(radar.data_smooth, 2));
widths = cell(1,size(radar.data_smooth, 2));
depths = cell(1, size(radar.data_smooth, 2));
depth_idx = cell(1, size(radar.data_smooth, 2));

for i = 1:size(radar.data_smooth, 2)
    data_i = radar.data_smooth(:,i);
    minProm = 0.05;                 % Prominence threshold for peaks
%     minProm = iqr(data_i);
%     minDist = max([0 min(dist_Ex)-3*dist_STD]);
    minDist = 0.12;                 % Min distance between peaks (in meters)
    
    % Find peaks in each trace based on requirements
    [~, peaks_idx_i, widths_i, Prom_i] = findpeaks(data_i, 'MinPeakProminence', ...
        minProm, 'MinPeakDistance', minDist/resolution);
%     peaks_i = zeros(length(data_i), 1);
%     peaks_i(peaks_idx_i) = Prom_i;
%     peaks(:,i) = peaks_i;

    % Add peak prominence and width values to relevent matrices
    peaks(peaks_idx_i,i) = Prom_i;
    peak_width(peaks_idx_i,i) = widths_i;
    
    % Add values to relevent cell arrays
    Proms{i} = Prom_i;
    widths{i} = widths_i;
    depths{i} = radar.depth(peaks_idx_i);
    depth_idx{i} = peaks_idx_i;
end

%%

% Define size of ~quasi bin confidence interval
err_bin = 6;
err_bin = minDist/resolution;

% Preallocate cell array for layer numbers and initialize values by
% assigning unique layer numbers to each peak in the first trace
Groups = Proms;
Groups{1} = uint32((1:length(Groups{1})))';

% Set value for next unique layer number
new_group = Groups{1}(end) + 1;

% Preallocate matrix for layer numbers and add initialized values for the
% first trace (col=1)
peak_group = zeros(size(peaks));
peak_group(depth_idx{1},1) = Groups{1};

for i = 2:size(peaks, 2)
    
    % Determine column bounds for the local search window for i based on 
    % the bin error window
    col_idx = [max([i-round(0.5*err_bin) 1]) i-1];
    
    for j = 1:length(Proms{i})
        
        % Determine cell array index for j in ith trace
        j_idx = depth_idx{i}(j);
        
        % Determine index values for the row boundaries of the local search
        % window of peak (i,j) based on the bin error window size and the
        % half-width of peak (i,j)
        row_idx = [max([j_idx-round(0.5*(err_bin+widths{i}(j))) 1]) ...
            min([j_idx+round(0.5*(err_bin+widths{i}(j))) size(peaks, 1)])];
        
        % Define local window to search for matching layer numbers
        peaks_local = peaks(row_idx(1):row_idx(2),col_idx(1):col_idx(2));
        
        % Find the row, col, and index values for peaks within the local
        % search window
        [local_row, local_col] = find(peaks_local);
        local_idx = sub2ind(size(peaks_local), local_row, local_col);
        
        % Calculate distances between peak (i,j) and peaks within the local
        % search window based on differences in depth, lateral distance,
        % and peak prominence magnitudes
        w_dist = sqrt((j_idx - (row_idx(1)+local_row-1)).^2 + ...
            (i - (col_idx(1)+local_col-1)).^2 + ...
            (Proms{i}(j) - peaks_local(local_idx)).^2);
        
        % Select the nearest neighbor to peak (i,j)
        [dist_min, dist_idx] = min(w_dist);
        
        % (I may add a tolerance in the future)
        if ~isempty(dist_idx)     % if j_min <= bin_res/2
            
            % If present, assign the neighest neighbor layer number to the 
            % peak (i,j) layer number in both cell array and matrix
            group_j = uint32(peak_group(row_idx(1)+local_row(dist_idx)-1,...
                col_idx(1)+local_col(dist_idx)-1));
            Groups{i}(j) = group_j;
            peak_group(j_idx,i) = group_j;
        else
            
            % If peak (i,j) does not have a nearest neighbor, assign a new
            % unqiue layer number to peak (i,j)
            Groups{i}(j) = new_group;
            peak_group(j_idx,i) = new_group;
            new_group = new_group + 1;
        end
    end
end

layers_idx = cell(1,new_group-1);
layers_val = zeros(1, length(layers_idx));
layer_peaks = zeros(size(peaks));
for i = 1:length(layers_idx)
    layers_idx{i} = find(peak_group==i);
    layers_val(i) = sum(peaks(layers_idx{i}));
    layer_peaks(layers_idx{i}) = layers_val(i);
end


age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);

P_50 = mean(mean(peaks(peaks>0)))*1000/mode(diff(radar.dist));
Po = 0.05;
K = 1;
r = log((K*Po/0.50-Po)/(K-Po))/-P_50;

ages = zeros([size(radar.data_smooth) Ndraw]);

for i = 1:size(layer_peaks, 2)
    
    peaks_i = layer_peaks(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar.depth(peaks_idx);
    
    % Probability of peak representing a year based on a logistic function
    % with rate (r) calculated above
    likelihood = K*Po./(Po + (K-Po)*exp(-r*peaks_i));
    
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end
    
    for j = 1:Ndraw
        depths_j = [0; depths_i(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        ages(:,i,j) = interp1(depths_j, yrs_j, radar.depth, 'linear', 'extrap');
    end
end

radar.age = ages;

%% Diagnostic figures

% Select random radar trace for comparison plots
i = randi(size(radar.data_smooth, 2));

% Plot radargram
figure('Position', [200 200 1500 800])
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
colorbar
xlabel('Distance along profile (m)')
ylabel('Depth (m)')
hold on
plot([radar.dist(i) radar.dist(i)], [0 radar.depth(end)], 'r', 'LineWidth', 2)
xlim([0 radar.dist(end)])
ylim([0 radar.depth(end)])
set(gca, 'Ydir', 'reverse', 'FontSize', 18)
hold off

% Compare radar signal at random trace to annual horizion selection at that
% trace
% layer_idx = logical([diff(floor(squeeze(radar.age(:,i,100)))); 0]);
layer_idx = logical([diff(floor(squeeze(radar.age(:,i,randi(Ndraw))))); 0]);
figure('Position', [50 50 500 1200])
hold on
plot(radar.data_smooth(:,i), radar.depth, 'r', 'LineWidth', 1.5)
scatter(radar.data_smooth(layer_idx,i), radar.depth(layer_idx), 'b<', 'filled')
xlim([-1.25 2])
xlabel('Radar Z-score')
ylabel('Depth (m)')
set(gca, 'Ydir', 'reverse')
hold off

% Compare radar signal at random trace to the mean annual horizions at that
% same random trace
layer_idx = logical([diff(floor(mean(squeeze(radar.age(:,i,:)), 2))); 0]);
figure('Position', [50 50 1200 500])
hold on
plot(radar.depth, radar.data_smooth(:,i))
plot([radar.depth(layer_idx) radar.depth(layer_idx)], ...
    [-1.5 1.5], 'k');
ylim([-1.5 1.5])


age_mean = mean(squeeze(radar.age(:,i,:)), 2);
age_ERR = 2*std(squeeze(radar.age(:,i,:)), [], 2);

figure
hold on
h1 = plot(core.depth, core.age, 'b', 'LineWidth', 2);
h2 = plot(radar.depth, age_mean, 'r', 'LineWidth', 2);
plot(radar.depth, age_mean + age_ERR, 'r--', 'LineWidth', 0.5)
plot(radar.depth, age_mean - age_ERR, 'r--', 'LineWidth', 0.5)
ylabel('Calendar Year')
xlabel('Depth (m)')
ylim([min([min(core.age) min(age_mean-age_ERR)]) max([max(core.age) max(age_mean)])])
legend([h1 h2], 'Core age (manual)', 'Radar age (automated)', 'Location', 'ne')
set(gca, 'FontSize', 10)
hold off

toc
