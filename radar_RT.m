function [radar] = radar_RT(radar_file, cores)

radar_type = isstruct(radar_file);

switch radar_type
    case false
        % Conversion to depth
        [radar] = radar_depth(radar_file, cores);
        
    case true
        % Conversion to depth
        [radar] = OIB_depth(radar_file, cores);
end

% % Determine if there are data break points (large gaps in data that would
% % necessitate data processing over segments of the whole data)
% distance = [0 diff(radar.dist)];
% dist_idx = distance >= 500;
% data_col = 1:size(radar.data_out, 2);
% data_endpts = [1 data_col(dist_idx)-1 length(data_col)];

% Find the mean response with depth in the radar data attributes across a
% given horizontal resolution (in meters)
horz_res = 25;
[radar] = radar_stack(radar, horz_res);

% Stationarize the radar response by differencing traces with a smoothing 
% spline
s = zeros(size(radar.data_stack));
for i = 1:size(s, 2)
    s(:,i) = csaps(radar.depth(:,i), radar.data_stack(:,i), 0.95, radar.depth(:,i));
end
radar_stat = radar.data_stack - s;

% Remove linear trend in variance (attentuation with depth) and convert to
% standardized values (z-score statistics)
radar_Z = zeros(size(radar_stat));
for i = 1:size(radar_stat, 2)
    data = radar_stat(:,i);
    % Frame length to define local variance
    half_frame = round(0.5*length(data)/5); 
    var0 = movvar(data, 2*half_frame, 'EndPoints', 'discard');
    x = (half_frame:length(data)-half_frame)';
    % Linear trend in variance
    EQ = polyfit(x, var0, 1);
    x_mod = (1:length(data))';
    mod = polyval(EQ, x_mod);
    % Standardize variance-corrected data
    radar_Z(:,i) = data./sqrt(abs(mod));
end

% Define the vertical resolution of the core data
core_res = 0.02;

% Define the cutoff depth for radar traces and find index of crossover
% depth
cutoff = 25;
depth_bott = floor(min([min(radar.depth(end,:)) cutoff]));

% Trim radar traces to cutoff depth and interpolate data to vertical scale
% of the firn cores
radarZ_interp = zeros(depth_bott/core_res+1, size(radar.data_stack, 2));
for i = 1:size(radar.data_stack, 2)
    depth_interp = (0:core_res:radar.depth(end,i));
    radarZ_i = interp1(radar.depth(:,i), radar_Z(:,i), depth_interp, 'pchip');
    radarZ_interp(:,i) = radarZ_i(1:size(radarZ_interp, 1));
end

% If manual layer picks are present, transform them to same depth and
% vertical scale as the interpolated radar data
if isfield(radar, 'man_layers')
    man_interp = zeros(size(radarZ_interp));
    for i = 1:size(radar.data_stack, 2)
%         depth_interp = (0:core_res:radar.depth(end,i))';
        layer_idx = logical(radar.man_layers(:,i));
        layer_num = radar.man_layers(layer_idx,i);
        man_depth_i = radar.depth(layer_idx,i);
        depth_idx = round(man_depth_i/core_res) + 1;
        man_interp(depth_idx(man_depth_i<=cutoff),i) = layer_num(man_depth_i<=cutoff);
    end
    radar.man_layers = man_interp;
end

% Assign structure output depth to interpolated depths
radar.depth = (0:core_res:depth_bott)';

% Smooth the laterally averaged radar traces with depth based on a 3rd
% order Savitzky-Golay filter with a window of 9 frames (~20 m)
radar.data_smooth = sgolayfilt(radarZ_interp, 3, 9);

% Clear unnecessary variables
clearvars -except file cores Ndraw radar horz_res core_res

%% Iterative radon transforms

depth_interval = 5;
dist_interval = 150;
depth_sz = round(0.5*depth_interval/core_res);
dist_sz = round(0.5*dist_interval/horz_res);

s_matrix = nan(size(radar.data_smooth));
ii = dist_sz+1:3:size(radar.data_smooth,2)-dist_sz;
jj = depth_sz+1:5:size(radar.data_smooth,1)-depth_sz;

for i = 1:length(ii)
    for j = 1:length(jj)
        data = radar.data_smooth(jj(j)-depth_sz:jj(j)+depth_sz,ii(i)-dist_sz:ii(i)+dist_sz);
        theta = 0:179;
        [R,~] = radon(data, theta);
        
        [R_max,idx] = max(R, [], 2);
        R_max = sum(R,2);
        R_max(R_max<0) = 0;
        weights = R_max/(sum(R_max));
        theta_max = sum(weights'.*theta(idx));
        s_matrix(jj(j),ii(i)) = -1*tand(theta_max-90);  
    end
    s_matrix(1,ii(i)) = 0;
    p = robustfit(1:length(radar.depth), s_matrix(:,ii(i)));  % Relation may be nonlinear
    s_matrix(end,ii(i)) = length(radar.depth)*p(2);
end


s_matrix(:,1) = s_matrix(:,ii(1));
s_matrix(:,end) = s_matrix(:,ii(end));


x = [1 ii size(radar.data_smooth,2)];
y = [1 jj size(radar.data_smooth,1)];
[X,Y] = meshgrid(radar.dist(x), radar.depth(y));
ss = s_matrix(y,x);
[Vx, Vy] = meshgrid(radar.dist, radar.depth);
Vq = interp2(X, Y, ss, Vx, Vy);
grad_smooth = imgaussfilt(Vq, 5);


%% Other tests with linear regression




% Diagnostic plot
ystart = 10:25:size(grad_smooth,1);
xstart = ones(1, length(ystart));
XY_raw = stream2(ones(size(grad_smooth)), grad_smooth, xstart, ystart, 1);
XY = XY_raw;
for k = 1:length(XY)
    XY{k}(:,1) = XY_raw{k}(:,1)*mean(diff(radar.dist));
    XY{k}(:,2) = XY_raw{k}(:,2)*.02;
end
figure
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
hold on
hlines = streamline(XY);
set(hlines, 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '--')
hold off



%%

end