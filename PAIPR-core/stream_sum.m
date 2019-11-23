function [stream_val, XY] = stream_sum(radar, IM_grad)
%STREAM_SUM Summary of this function goes here
%   Detailed explanation goes here

ystart = 1:size(IM_grad,1);
xstart = repmat(round(length(radar.dist)/2), 1, length(ystart));
% xstart = ones(1, length(ystart));
XY_left = stream2(-1*ones(size(IM_grad)), -1*IM_grad, xstart, ystart, 1);
XY_right = stream2(ones(size(IM_grad)), IM_grad, xstart+1, ystart, 1);
XY_raw = cellfun(@(x,y) [flipud(x);y], XY_left, XY_right, ...
    'UniformOutput', false);

XY = cellfun(@(x) sub2ind(size(radar.data_smooth), round(x(:,2)), round(x(:,1))), ...
    XY_raw, 'UniformOutput', false);
horz_res = round(mean(diff(radar.dist)));

% XY_smooth = XY;
% for i = 1:length(XY_smooth)
%     for j = 1:length(XY{i})
%         
%         % Get subscripts
%         [r,c] = ind2sub(size(radar.data_smooth), XY{i}(j));
%         
%         % Define smoothing window box
%         x_bnd = [max([1 c-5]) min([size(radar.data_smooth,2) c+5])];
%         y_bnd = [max([1 r-5]) min([size(radar.data_smooth,1) r+5])];
%         j_box = radar.data_smooth(y_bnd(1):y_bnd(2),x_bnd(1):x_bnd(2));
%         
%         XY_smooth{i}(j) = mean(j_box(:));
%     end
% end

stream_val = horz_res*cellfun(@(x) sum(radar.data_smooth(x)), XY);


% % Diagnostic plot
% figure
% imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
% hold on
% XY_plot = XY_raw;
% for k = 1:length(XY)
%     XY_plot{k}(:,1) = XY_raw{k}(:,1)*mean(diff(radar.dist));
%     XY_plot{k}(:,2) = XY_raw{k}(:,2)*.02;
% end
% hlines = streamline(XY_plot(1:10:length(XY_plot)));
% set(hlines, 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '--')
% hold off

end

