function [stream_val, XY] = stream_sum(radar, IM_grad, xstart_idx)
%STREAM_SUM Summary of this function goes here
%   Detailed explanation goes here

ystart_vec = 1:size(IM_grad,1);


% XY_sets = cell(1,length(xstart_idx));
% for i=1:length(xstart_idx)
%     
%     xstart = repmat(xstart_idx(i), 1, length(ystart_vec));
%     XY_left = stream2(-1*ones(size(IM_grad)), -1*IM_grad, ...
%         xstart, ystart_vec, 1);
%     XY_right = stream2(ones(size(IM_grad)), IM_grad, ...
%         xstart+1, ystart_vec, 1);
%     XY = cellfun(@(x,y) [flipud(x);y], XY_left, XY_right, ...
%         'UniformOutput', false);
%     
%     XY_sets{i} = cellfun(@(x) sub2ind(size(radar.data_smooth), ...
%         round(x(:,2)), round(x(:,1))), XY, 'UniformOutput', false);
%     
% end
% 
% XY = cell(sum(cellfun(@length, XY_sets)),1);
% XY(1:length(XY_sets{1})) = XY_sets{1};
% XY_loop = XY_sets(2:end);
% node_idx = length(XY_sets{1})+1;
% for i=1:length(XY_loop)
%     
%     XY_i = XY_loop{i}';
%     for j=1:length(XY_i)
%         XY_clip = XY(~cellfun(@isempty,XY));
%         j_intersects = zeros(length(XY_clip),1);
%         for k = 1:length(j_intersects)
%             j_intersects(k) = length(intersect(XY_i{j},XY_clip{k}));
%         end
%         [j_max,j_idx] = max(j_intersects);
%         if j_max < 1
%             XY{node_idx} = XY_i{j};
%             node_idx = node_idx+1;
%         end
%     end
% end
   

xstart = repmat(xstart_idx, 1, length(ystart_vec));
XY_left = stream2(-1*ones(size(IM_grad)), -1*IM_grad, ...
    xstart, ystart_vec, 1);
XY_right = stream2(ones(size(IM_grad)), IM_grad, ...
    xstart+1, ystart_vec, 1);
XY_raw = cellfun(@(x,y) [flipud(x);y], XY_left, XY_right, ...
    'UniformOutput', false);

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

XY = cell(1,length(XY_raw));
for i=1:length(XY)
    
    col = unique(round(XY_raw{i}(:,1)));
    col(col<1) = 1;
    col(col>size(radar.data_smooth,2)) = size(radar.data_smooth,2);
    row = round(interp1(XY_raw{i}(:,1),XY_raw{i}(:,2), col, ...
        'linear', 'extrap'));
    row(row<1) = 1;
    row(row>size(radar.data_smooth,1)) = size(radar.data_smooth,1);
    XY{i} = sub2ind(size(radar.data_smooth), row, col);
end

% XY = cellfun(@(x) sub2ind(size(radar.data_smooth), ...
%     round(x(:,2)), round(x(:,1))), XY_raw, 'UniformOutput', false);

horz_res = round(mean(diff(radar.dist)));
stream_val = horz_res*cellfun(@(x) sum(radar.data_smooth(x)), XY);

end

