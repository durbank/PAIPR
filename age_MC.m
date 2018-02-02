function [radar, core_syn] = age_MC(file, cores, Ndraw)

% Conversion to depth
[radar, core_syn] = radar_depth(file, cores);

% Stationarize the radar response using a smoothing spline
s = zeros(size(radar.data_out));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth, radar.data_out(:,i), 0.95, radar.depth);
end
radar_stat = radar.data_out - s;

% Remove linear trend in variance (attentuation with depth) and convert to
% z-score statistics
radar_Z = zeros(size(radar_stat));
for i = 1:size(radar_stat, 2)
    data = radar_stat(:,i);
    half_frame = round(0.5*length(data)/5);
    var0 = movvar(data, 2*half_frame, 'EndPoints', 'discard');
    x = (half_frame:length(data)-half_frame)';
    EQ = polyfit(x, var0, 1);
    x_mod = (1:length(data))';
    mod = EQ(1)*x_mod+EQ(2);
    radar_Z(:,i) = zscore(data./sqrt(abs(mod)));
end

% % Scale radar data by converting to power
% radar_power = 20*log10(radar.data_out);
% 
% % Stationarize the radar power using a smoothing spline
% s = zeros(size(radar_power));
% for i = 1:size(s, 2)
%     s(:,i) = csaps(radar.depth, radar_power(:,i), 0.95, radar.depth);
% end
% radar_stat = radar_power - s;
% 
% % Assign radar_Z as radar_stat (for ease of switching to-from log-power
% radar_Z = radar_stat;

% EX = randi(size(radar.data_out, 2));
% % Diagnostic figure for stationarity
% figure
% hold on
% plot(radar.data_out(:,EX), radar.depth)
% plot(s(:,EX), radar.depth)
% xlabel('Radar return strength')
% ylabel('Depth (m)')
% legend('Original data', 'Stationarizing spline', 'location', 'nw')
% set(gca, 'Ydir', 'reverse')
% hold off

% Resample and interpolate radar depth and data to a constant depth
% resolution using a shape-preserving cubic interpolation
resolution = 0.02;       % Depth resolution in meters
radar.depth_interp = (0:resolution:radar.depth(end))';
radarZ_resamp = interp1(radar.depth, radar_Z, radar.depth_interp, 'pchip');

% % Diagnostic figure
% figure
% plot(radarZ_resamp(:,EX), radar.depth_interp, 'k')
% xlabel('Radar return strength')
% ylabel('Depth (m)')
% set(gca, 'Ydir', 'reverse')

% Find the mean response with depth in the resampled radar data across a
% given lateral distance 'window' (lateral bin size, where 1 bin ~ 2.5 m)
window = 50;
radar_mean = movmean(radarZ_resamp, window, 2);

% % Diagnostic figure
% figure
% plot(radar_mean(:,EX), radar.depth_interp, 'm')
% xlabel('Radar return strength')
% ylabel('Depth (m)')
% set(gca, 'Ydir', 'reverse')

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~20 m)
radar.data_smooth = sgolayfilt(radar_mean, 3, 9);

% % Diagnostic figure
% figure
% plot(radar.data_smooth(:,EX), radar.depth_interp, 'r')
% xlabel('Radar return strength')
% ylabel('Depth (m)')
% set(gca, 'Ydir', 'reverse')

% Year associated with the first pick of the algorithm
% age_top = round(radar.collect_date);
age_top = radar.collect_date;
yr_pick1 = ceil(age_top-1);

% Indices of annual horizions in core
yr_idx = logical([diff(floor(core_syn.age)); 0]);

% Difference in depth between annual horizions in core
depth_diff = diff(core_syn.depth(yr_idx));
depth_btw_yr = [depth_diff; mean(depth_diff)];

% Model of distance between annual horizions with depth in core
f = polyfit(core_syn.depth(yr_idx), depth_btw_yr, 1);

% Expected thickness (and standard deviation) of annual layers 
% throughout the depth of the core
dist_Ex = f(1)*radar.depth_interp + f(2);
dist_STD = std(depth_btw_yr);

% Estimate the depth of annual peaks within each smoothed radar trace
ages = zeros([size(radar.data_smooth) Ndraw]);
data = radar.data_smooth;
depth_i = radar.depth_interp;
% Prom_all = cell(1, size(data, 2));
parfor i = 1:size(radar.data_smooth, 2)
    data_i = data(:,i);
    minProm = 0.01;
    minDist = abs(min(dist_Ex) - 3*dist_STD);
    [~, depths_peaks, ~, Prom] = findpeaks(data_i, depth_i, ...
        'MinPeakProminence', minProm, 'MinPeakDistance', minDist);
    
    % Remove the surface horizon peak from the data
    [~,max_idx] = max(Prom);
    Prom(max_idx) = [];
    depths_peaks(max_idx) = [];
    
    % Solve for r so that prominence with 0.50 probability is equal to 
    % interquartile range of the smoothed data (may need to be optimized to
    % observations using a scaling constant)
    Po = 0.01;
    K = 1;
    P_50 = 1*iqr(data_i);
%     P_50 = median(Prom);
    r = log((K*Po/0.50-Po)/(K-Po))/-P_50;
    
    % Probability of peak representing a year based on a logistic function
    % with rate (r) calculated above
    P = K*Po./(Po + (K-Po)*exp(-r*Prom));
%     plot(Prom, P, 'bo')

%%% EXPERIMENTAL %%%
    % Calculate based on fitted Gamma distribution and corresponding
    % cummulative density function
%     idx = Prom>quantile(Prom, 0.25);
%     Prom = Prom(idx);
%     depths_i = depths_i(idx);
%     G_params = gamfit(Prom);
%     G_pdf = pdf('Gamma', 0:0.01:ceil(max(Prom)), G_params(1), G_params(2));
%     P = cdf('Gamma', Prom, G_params(1), G_params(2));
%%% EXPERIMENTAL %%%

    % Calculate distance weights based off of the distances to the nearest
    % four peaks (in future should scale the effect based off the variance 
    % in the expected distance, dist_Ex_i)
    [~, Prom_dist] = knnsearch(depths_peaks, depths_peaks, 'K', 5);
    dist_Ex_i = f(1)*depths_peaks + f(2);
%     w = mean(Prom_dist(:,2:end)-dist_Ex_i, 2)./dist_Ex_i;
    w = sum([0.5 0.5 0.25 0.25].*(Prom_dist(:,2:end)-dist_Ex_i), 2)./dist_Ex_i;
    Pw = P.^(1-w);
    P = Pw;
    
    yr_idx = zeros(length(depths_peaks), Ndraw);
    for j = 1:length(depths_peaks)
        R = rand(Ndraw, 1) <= P(j);
        yr_idx(j,:) = R;
    end
    
    age_i = zeros(length(data_i), Ndraw);
    for j = 1:Ndraw
        depths_j = [0; depths_peaks(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
%         depths_j = depths_peaks(logical(yr_idx(:,j)));
%         yrs_j = (age_top:-1:age_top-length(depths_j)+1)';
        age_i(:,j) = interp1(depths_j, yrs_j, depth_i, 'linear', 'extrap');
    end
    ages(:,i,:) = age_i;
end
radar.age = ages;

radar = rmfield(radar, {'time_trace', 'arr_layers', 'TWTT'});
end


