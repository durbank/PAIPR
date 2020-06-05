% Function to evaluate echogram QC image to determine if there are issues
% regarding the quality of the echogram data

function [med_val, threshold_val, QC_flag, depth_idx, yr_cutoff] = ...
    QC_check(mdata, warning_lim, thresh_perc)

QC_image = mdata.IM_QC;


% % figure
% % imagesc(radar.data_smooth, [-3 3])
% % figure
% % plot(radar.depth, mean(QC_image, 2))

% figure
% imagesc(QC_image, [0 1])
% colorbar

% Get array of all nonzero QC values
score_vals = QC_image(QC_image > 0);

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

% Get array of median ages
ages = median(mdata.ages, 3);

% Check if QC passes both overall limit and limit throughout depth
if (threshold_val <= warning_lim) && (max(depth_vals) <= warning_lim)
    
    % Set QC flag
    QC_flag = 0;
    
    % Set depth cutoff index for "good" data to bottom of image
    depth_idx = size(QC_image,1);

    % Set year cutoff to floor of oldest age
    yr_cutoff = floor(min(ages(depth_idx,:)));
    
else
    
    % Find depth cutoff index for "good" data
    depth_idx = find(depth_vals > warning_lim, 1, 'first');
    if isempty(depth_idx)
        % Set depth cutoff index for "good" data to bottom of image
        depth_idx = size(QC_image,1);
        
        % Expand cutoffs for warning limit and % threshold by 20%
        warning_lim = warning_lim + 0.20*warning_lim;
        thresh_perc = thresh_perc + 0.20*thresh_perc;
    end
    
    % QC image above depth cutoff
    im_cut = QC_image(1:depth_idx, :);
    % Get array of all nonzero QC values
    vals_cut = im_cut(im_cut > 0);
    % Generate cummulative distribution for QC values
    x_cdf_cut = (1:length(vals_cut))/length(vals_cut);
    score_cdf_cut = sort(vals_cut);
    % Determine location where cdf exceeds cutoff threshold percentile
    cut_idx_cut = find(x_cdf_cut>(1-thresh_perc), 1, 'first');
    % Determine cdf value at cutoff threshold percentile
    thresh_val_cut = score_cdf_cut(cut_idx_cut);
    
    
    % Determine surface age and age of last "good" year
    yr_top = floor(max(ages(1,:)));
    yr_cutoff = ceil(mean(ages(depth_idx,:)));
    
    % Determine depth index for 2 m below surface
    depth_2m = round(2/mode(diff(mdata.depth)));
    
    % Check if <10% of data bad above depth_idx, contains at least 1 year
    % of reliable data, and reliable data extends at least 2m below surface
    if (thresh_val_cut <= warning_lim) && ...
            (yr_cutoff < yr_top) && (depth_idx > depth_2m)
        
        % Assign flag of data reliable to certain depth
        QC_flag = 1;
        
    else
        
        % Assign flag as all data questionable
        QC_flag = 2;
    end
    
end


% % Flag image if either % threshold value exceeds warning limit or if mean
% % QC value exceeds warning limit beyond a certain depth
% if (threshold_val >= warning_lim) || (max(depth_vals) >= warning_lim)
%     QC_flag = true;
%     
%     % Find depth cutoff index for "good" data
%     depth_idx = find(depth_vals > warning_lim, 1, 'first');
%     if isempty(depth_idx)
%         depth_idx = size(QC_image,1);
%     end
%     
% else
%     QC_flag = false;
%     
%     % Set depth cutoff index for "good" data to bottom of image
%     depth_idx = size(QC_image,1);
% end
% yr_cutoff = ceil(mean(ages(depth_idx,:)));


end