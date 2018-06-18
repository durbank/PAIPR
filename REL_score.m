% Function that looks at the similiarity in relative row position between
% different picked layers (a measure of how effectively/consistently the
% method picked layers in the given data)

function [RMSE_globe, RSE_layer, depth_slope] = REL_score(peaks, layers_idx)

[~, global_idx] = sort(cellfun(@length, layers_idx), 'descend');
layers_global = {sub2ind(size(peaks), ones(size(peaks,2),1), ...
    (1:size(peaks,2))') layers_idx{global_idx(1:24)}};

[r_all, c_all] = cellfun(@(x) ind2sub(size(peaks), x), layers_global, 'UniformOutput', false);

col_idx = [max(cellfun(@min, c_all)) min(cellfun(@max, c_all))];
cols = col_idx(1):col_idx(2);

% rows = cellfun(@(x) x(col_idx(1):col_idx(2)), r_all, 'UniformOutput', false);
rows = cell(1, length(c_all));
for i = 1:length(rows)
    r_idx = c_all{i} >= col_idx(1) & c_all{i} <= col_idx(2);
    rows{i} = r_all{i}(r_idx);
end



[row_depth, sort_idx] = sort(cellfun(@mean, rows'));
row_mat = (cell2mat(cellfun(@(x) movmean(x - mean(x), 10), rows, 'UniformOutput', false)))';
row_dev = row_mat(sort_idx,:);

depth_slope = zeros(1, length(cols));
RMSE_globe = zeros(1, length(cols));
for i = 1:length(cols)
%     [p] = polyfit(row_depth, row_dev(:,i), 1);
%     depth_slope(i) = p(1);
%     RMSE(i) = sqrt(mean((row_dev(:,i) - polyval(p, row_depth)).^2));
    [p, stats] = robustfit(row_depth, row_dev(:,i));
    depth_slope(i) = p(2);
    RMSE_globe(i) = stats.robust_s;
end


mdl = (depth_slope+1).*row_depth;

figure
hold on
for i = 1:length(rows)
plot(cols, rows{i}, 'LineWidth', 2)
end
set(gca, 'Ydir', 'reverse')
plot(cols, mdl, 'k')
hold off

% m = zeros(1, length(cols));
% RMSE = zeros(1, length(cols));
% for i = 1:length(cols)
%     mdl = fitlm(row_depth, row_dev(:,i));
%     m(i) = mdl.Coefficients.Estimate(2);
%     RMSE(i) = mdl.RMSE;
% end

% p = row_depth\row_dev;
% row_diff = mean(rows{sort_idx(15)} - (p*row_depth(15)+row_depth(15))');


%% Individual reliability scores

[r_all, c_all] = cellfun(@(x) ind2sub(size(peaks), x), layers_idx, 'UniformOutput', false);
row_depth = cellfun(@mean, r_all');

RSE_layer = cell(1, length(layers_idx));
for i = 1:length(RSE_layer)
    mdl = (depth_slope(c_all{i}(1):c_all{i}(end))+1).*row_depth(i);
    RSE_layer{i} = sqrt((r_all{i} - mdl').^2);
end





end