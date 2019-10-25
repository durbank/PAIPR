% Function to perform local radon transforms on an input echogram image. 
% Output is a matrix (same dimensions as input image) with the estimated
% layer slope at each bin location in the image.

function [IM_gradients] = radar_gradient(radar, vert_res, horz_res)

% Define depth/distance intervals over which to perform radon transforms
% (in meters), and calculate data matrix window size (in data bins)
depth_interval = 4;
dist_interval = 250;
depth_sz = round(0.5*depth_interval/vert_res);
dist_sz = round(0.5*dist_interval/horz_res);

% Define (overlapping) window center points for iterative radon transforms
ii = dist_sz+1:4:size(radar.data_smooth,2)-dist_sz;
jj = depth_sz+1:round(depth_sz/10):size(radar.data_smooth,1)-depth_sz;

% Preallocate matrix for estimated layer slopes
s_matrix = nan(size(radar.data_smooth));

% Iteratively calculate local layer gradients (data array units) for
% overlapping local windows within the echogram
for i = 1:length(ii)
    for j = 1:length(jj)
        
        % Define the local data window
        data_i = radar.data_smooth(jj(j)-depth_sz:...
            jj(j)+depth_sz,ii(i)-dist_sz:ii(i)+dist_sz);
        
        % Angles over which to perform radon transform
        theta = 0:179;
        
        % Calculate transform
        [R,~] = radon(data_i, theta);
        
        % Find the sum of all positive transformed values for each angle
        R_sum = sum(R.*(R>0));
        
        % Find the maximum positive sum value, and assign it's
        % corresponding angle as the dominant angle of layer gradients
        [~,max_idx] = max(R_sum);
        theta_max = theta(max_idx);

        % Transform dominant angle to true layer gradient and assign to
        % preallocated matrix
        s_matrix(jj(j),ii(i)) = -1*tand(theta_max-90);  
    end
    
    % Define the surface layer to have a gradient of 0
    s_matrix(1,ii(i)) = 0;
    
    % Extrapolate the layer gradients at the base of the echogram based on
    % the robust linear regression of the slopes of overlying layers
    p = robustfit(1:size(radar.depth,1), s_matrix(:,ii(i)));
    s_matrix(end,ii(i)) = size(radar.depth,1)*p(2);
end

% Extrapolate gradient values for echogram edges based on nearest non-NaN
% values
s_matrix(:,1) = s_matrix(:,ii(1));
s_matrix(:,end) = s_matrix(:,ii(end));

% Update the along-trace and depth locations of non-NaN values in slope
% matrix
x = [1 ii size(radar.data_smooth,2)];
y = [1 jj size(radar.data_smooth,1)];

% Interpolate layer gradients for every trace and depth point
[X,Y] = meshgrid(radar.dist(x), radar.depth(y));
ss = s_matrix(y,x);
[Vx, Vy] = meshgrid(radar.dist, radar.depth);
IM_gradients = interp2(X, Y, ss, Vx, Vy);


% % Diagnostic plot
% ystart = 1:25:size(IM_gradients,1);
% xstart = ones(1, length(ystart));
% XY_raw = stream2(ones(size(IM_gradients)), IM_gradients, ...
%     xstart, ystart, 1);
% XY = XY_raw;
% for k = 1:length(XY)
%     XY{k}(:,1) = XY_raw{k}(:,1)*mean(diff(radar.dist));
%     XY{k}(:,2) = XY_raw{k}(:,2)*.02;
% end
% figure
% imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
% hold on
% hlines = streamline(XY);
% set(hlines, 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '--')
% hold off

end