% Function that looks at the similiarity in relative row position between
% different picked layers (a measure of how effectively/consistently the
% method picked layers in the given data)

function [RMSE_globe, RSE_layer, depth_slope] = REL_score2(peaks, layers, layers_idx)


[r_all, c_all] = cellfun(@(x) ind2sub(size(peaks), x), ...
    {sub2ind(size(peaks), ones(size(peaks,2),1), (1:size(peaks,2))') ...
    layers_idx}, 'UniformOutput', false);



seg_length = round(mean(cellfun(@length, layers_idx)));
col_idx = [1:seg_length:size(peaks,2) size(peaks,2)];

depth_slope = zeros(1, size(peaks,2));
depth_cells = cell(1, size(peaks,2));
RMSE_globe = zeros(1, size(peaks,2));

figure
hold on
set(gca, 'Ydir', 'reverse')

for i = 1:length(col_idx)-1
    
    cols = col_idx(i):col_idx(i+1);
    layers_set = {sub2ind(size(peaks), ones(length(cols),1), (layers_idx(cellfun(@(x) min(x) <= col_idx(i) & max(x) >= col_idx(i+1), c_all));
    [r_set, c_set] = cellfun(@(x) ind2sub(size(peaks), x), layers_set, 'UniformOutput', false);
    
    rows = cell(1, length(c_set));
    for j = 1:length(rows)
        r_idx = c_set{j} >= col_idx(i) & c_set{j} <= col_idx(i+1);
        rows{j} = r_set{j}(r_idx);
    end
    
    
    [row_depth, sort_idx] = sort(cellfun(@mean, rows'));
    row_mat = (cell2mat(cellfun(@(x) x - mean(x), rows, 'UniformOutput', false)))';
    row_dev = row_mat(sort_idx,:);
    
    
    for j = 1:length(cols)
        %     [p] = polyfit(row_depth, row_dev(:,i), 1);
        %     depth_slope(i) = p(1);
        %     RMSE(i) = sqrt(mean((row_dev(:,i) - polyval(p, row_depth)).^2));
        [p, stats] = robustfit(row_depth, row_dev(:,j));
        depth_slope(col_idx(i)+j-1) = p(2);
        RMSE_globe(col_idx(i)+j-1) = stats.robust_s;
    end
    
    depth_cells(col_idx(i):col_idx(i+1)) = {row_depth};
    
    
    for j = 1:length(rows)
        plot(cols, rows{j}, 'LineWidth', 2)
    end
    plot(cols, mdl, 'k')
    
end

hold off

    




end