function [radar] = calc_SWE(radar, rho_data, Ndraw)


rho_mod = zeros(size(radar.data_smooth));
rho_std = zeros(size(radar.data_smooth));
for i = 1:size(rho_mod,2)
    
    depth_i = rho_data.Data{i}.Depth;
    mod_i = rho_data.Data{i}.pred_mean;
    std_i = rho_data.Data{i}.pred_sd;
    rho_mod(:,i) = interp1(depth_i, mod_i, radar.depth);
    rho_std(:,i) = interp1(depth_i, std_i, radar.depth);
end

%% Calculate annual accumulation for each radar trace

% Calculate accumulation at each depth interval (0.02 m) with simulated
% noise based on the variance in core density (also converts to mm w.e.)
noise_rho = 1000*repmat(rho_std, 1, 1, Ndraw).*randn(size(radar.ages));
accum_dt = 0.02*(1000*repmat(rho_mod, 1, 1, Ndraw) + noise_rho);

% Linearly interpolate traces with missing or problematic age-depth scales
% sim_nan = all(isnan(radar.ages), 3);
% trace_nan = any(sim_nan);
radar.ages = fillmissing(radar.ages, 'linear', 2);

% Define age-depth profile based on the median of all MC age profiles and 
% the std dev with depth of all MC age profiles for each trace (avoids  
% integer year jumps in accumulation estimates and non-monotonically 
% decreasing ages)
% age_std = squeeze(std(radar.ages, [], 3));
% age_noise = randn(1, 1, Ndraw).*age_std;
% ages = repmat(median(radar.ages, 3), 1, 1, Ndraw) + age_noise;
ages = radar.ages;

% Preallocation of cell arrays for accumulation years and annual
% accumulation rate
accum_yr = cell(1, size(radar.ages, 2));
accum = cell(1, size(radar.ages, 2));
for i = 1:size(accum, 2)
    
    % Define calendar age of the top of the first year with complete
    % accumulation data and earliest whole year with observations within the
    % data set
    yr_top = floor(max(ages(1,i,:)));
    yr_end = ceil(min(ages(end,i,:)));
    
    % Define initial accumulation year vector (will be iteratively modified at
    % each trace in for loop)
    accum_yr_init = (yr_top-1:-1:yr_end)';
    
    % Preallocate vector for accumulation rate at ith trace
    accum_i = nan(length(accum_yr_init), Ndraw);
%     accum_i = zeros(length(accum_yr_init), Ndraw);
    for j = 1:Ndraw
        
        % Calculate indices of integer ages for jth simulation of the ith
        % trace (with added noise from uncertainty in exact point in time
        % of the accumulation peak, using a std dev of 1 month)
        years_j = ages(:,i,j);
        yr_idx = logical([diff(floor(years_j)); 0]);
        yr_loc = find(yr_idx);
        loc_temp = yr_loc;
        loc_temp(2:end-1) = yr_loc(2:end-1) + ...
            round(movmean(1*diff(yr_loc(1:end-1))/12, 2).*...
            randn(length(yr_loc)-2, 1));
        loc_idx = loc_temp<1;
        loc_temp(loc_idx) = yr_loc(loc_idx);
        yr_loc = sort(loc_temp, 'ascend');
        yr_loc(yr_loc>size(accum_dt, 1)) = size(accum_dt, 1);
        
        % Integrate accumulation at each depth point for each whole year in
        % the jth simulation of the ith trace
        n_length = length(yr_loc) - 1;
        accum_j = zeros(n_length, 1);
        for n = 1:n_length
            accum_j(n) = sum(accum_dt(yr_loc(n)+1:yr_loc(n+1),i,j));
        end
        accum_i(1:n_length,j) = accum_j;
    end
    
    % Output results to respective arrays
    perc_num = sum(~isnan(accum_i),2)./Ndraw;
    accum_idx = find(perc_num >= 0.67, 1, 'last');
%     accum_idx = find(all(accum_i, 2), 1, 'last');
    accum_clip = accum_i(1:accum_idx,:);
    accum_yr{i} = accum_yr_init(1:accum_idx);
    accum{i} = accum_clip;
end


%% Output results to radar and core structures

radar.SMB_yr = accum_yr;
radar.SMB = accum;

end