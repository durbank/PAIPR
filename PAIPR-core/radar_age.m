% Function to calculate layer likelihood scores and age-depth profile
% distributions (based on likelihood scores)

function [radar] = radar_age(radar, r, k, Ndraw)

% % Find horizontal resolution of input radar echogram
% horz_res = mean(diff(radar.dist));
% 
% % Calculate continuous layer distances for each layer (accounting for 
% % lateral size of stacked radar trace bins)
% layers_dist = cellfun(@(x) numel(x)*horz_res, radar.layers);
% 
% % Find the mean peak prominence (used to scale the prominence-distance
% % results)
% % peak_w = 1/mean(radar.peaks(radar.peaks>0));
% % dist_w = 1/(size(radar.data_smooth,2)*horz_res);
% peak_w = 1;
% dist_w = 1;
% 
% % Map layer prominence-distance values to the location within the radar
% % matrix of the ith layer
% layer_DB = zeros(size(radar.peaks));
% for i = 1:length(radar.layers)
%     layer_DB(radar.layers{i}) = peak_w*dist_w*...
%         radar.peaks(radar.layers{i}).*layers_dist(i);
% end
% 
% radar.DB = layer_DB;

%% Assign layer likelihood scores and estimate age-depth scales

% Preallocate arrays for layer likelihoods and anges
ages = zeros([size(radar.data_smooth) Ndraw]);
radar.likelihood = zeros(size(radar.data_smooth));
err_out = [];
for i = 1:size(radar.DB, 2)
    
    % Define surface age and the year associated with the first pick of the
    % algorithm
    yr_vec = datevec(radar.collect_time(i));
    yr_pick1 = yr_vec(1);
    age_top = yr_vec(1) + (30*yr_vec(2)+yr_vec(3))/365;
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = radar.DB(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar.depth(peaks_idx);
    
    likelihood = 1./(1+exp(r*peaks_i + k));
    radar.likelihood(peaks_idx,i) = likelihood;
    
    % Assign MC simulation annual layer presence based on layer likelihood
    % values
    yr_idx = zeros(length(depths_i), Ndraw);
    for j = 1:length(depths_i)
        R = rand(Ndraw, 1) <= likelihood(j);
        yr_idx(j,:) = R;
    end 

    for j = 1:Ndraw
        depths_j = [0; depths_i(logical(yr_idx(:,j)))];
        yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
        try
            ages(:,i,j) = interp1(depths_j, yrs_j, radar.depth, ...
                'linear', 'extrap');
        catch
            sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
            err_out = [err_out j];
        end
    end
    if ~isempty(err_out)
        ages(:,i,err_out) = repmat(sum(squeeze(ages(:,i,:)), 2)./...
            sum(squeeze(ages(:,i,:))~=0, 2), 1, length(err_out));
    end
    err_out = [];
end

radar.ages = ages;

end