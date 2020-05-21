% Function to evaluate echogram QC image to determine if there are issues
% regarding the quality of the echogram data

function [med_val, threshold_val, QC_flag, depth_idx] = QC_check(...
    QC_image, warning_lim, thresh_perc)

% % QC_image = radar.im_QC;
% % figure
% % imagesc(radar.data_smooth, [-3 3])
% % figure
% % plot(radar.depth, mean(QC_image, 2))



% figure
% imagesc(QC_image, [0 1])
% colorbar

% Get array of all nonzero QC values
score_vals = QC_image(QC_image >0);

% [f, xi] = ksdensity(score_vals);
% figure
% plot(xi, f)
% xlim([0 1])

% Find the median QC value for the image
med_val = nanmedian(score_vals);

% Generate cummulative distribution for QC values
x_cdf = (1:length(score_vals))/length(score_vals);
score_cdf = sort(score_vals);

% figure
% plot(x_cdf, score_cdf)

% Determine location where cdf exceeds cutoff threshold percentile
cut_idx = find(x_cdf>(1-thresh_perc), 1, 'first');

% Determine cdf value at cutoff threshold percentile
threshold_val = score_cdf(cut_idx);

% Find the mean QC value with depth in image
depth_vals = nanmean(QC_image, 2);

% Flag image if either % threshold value exceeds warning limit or if mean
% QC value exceeds warning limit beyond a certain depth
if (threshold_val >= warning_lim) || (max(depth_vals) >= warning_lim)
    QC_flag = true;
    
    % Find depth cutoff index for "good" data
    depth_idx = find(depth_vals > warning_lim, 1, 'first');
    if isempty(depth_idx)
        depth_idx = size(QC_image,1);
    end
    
else
    QC_flag = false;
    
    % Set depth cutoff index for "good" data to bottom of image
    depth_idx = size(QC_image,1);
end

end