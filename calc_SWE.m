function [radar] = calc_SWE(radar, Ndraw)

% Generate modelled std of density with depth for radar data
rho_std = sqrt((radar.rho_var(1,:)-radar.rho_var(2,:))./...
    (radar.rho_var(3,:).^radar.depth) + radar.rho_var(2,:));

% Generate density with depth model for radar data
rho_mod = radar.rho_coeff(1,:).*radar.depth.^radar.rho_coeff(2,:) + ...
    radar.rho_coeff(3,:);

%% Calculate annual accumulation for each radar trace

% % Calculate accumulation at each depth interval (0.02 m) with simulated
% % noise based on the variance in core density
noise_rho = 1000*repmat(rho_std, 1, 1, Ndraw).*randn(size(radar.ages));
accum_dt = 0.02*(1000*repmat(rho_mod, 1, 1, Ndraw) + noise_rho);

% Define age-depth profile as the mean of all MC age profiles for each
% trace (avoids integer year jumps in accumulation estimates and
% non-monotonically decreasing ages)
age_std = squeeze(std(radar.ages, [], 3));
age_noise = randn(1, 1, Ndraw).*age_std;
ages = repmat(median(radar.ages, 3), 1, 1, Ndraw) + age_noise;
% ages = radar.age;



% Define first year with complete accumulation data and earliest year with
% observations within the data set
yr_top = floor(max(max(ages(1,:,:))));
yr_end = ceil(min(min(ages(end,:,:))));

% Define initial accumulation year vector (will be iteratively modified at
% each trace in for loop)
accum_yr_init = (yr_top-1:-1:yr_end)';

% Preallocation of cell arrays for accumulation years and annual
% accumulation rate
accum_yr = cell(1, size(radar.ages, 2));
accum = cell(1, size(radar.ages, 2));
for i = 1:size(accum, 2)
    
    % Preallocate vector for accumulation rate at ith trace
    accum_i = zeros(length(accum_yr_init), Ndraw);
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
        yr_loc = loc_temp;
        
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
    accum_idx = find(all(accum_i, 2), 1, 'last');
    accum_clip = accum_i(1:accum_idx,:);
    accum_yr{i} = accum_yr_init(1:accum_idx);
    accum{i} = accum_clip;
end


%% Output results to radar and core structures

radar.SMB_yr = accum_yr;
radar.SMB = accum;

end