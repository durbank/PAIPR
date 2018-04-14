function [radar] = calc_SWE(radar, Ndraw)

% Define first year with complete accumulation data and earliest year with
% observations within the data set
yr_top = floor(max(max(radar.age(2,:,:))));
yr_end = ceil(min(min(radar.age(end,:,:))));

% Generate modelled std of density with depth for radar data
rho_std = sqrt((radar.rho_var(1,:)-radar.rho_var(2,:))./...
    (radar.rho_var(3,:).^radar.depth) + radar.rho_var(2,:));

% Generate density with depth model for radar data
rho_mod = radar.rho_coeff(1,:).*radar.depth.^radar.rho_coeff(2,:) + ...
    radar.rho_coeff(3,:);

%% Calculate annual accumulation core each radar trace

% % Calculate accumulation at each depth interval (0.02 m) with simulated
% % noise based on the variance in core density
noise_rho = 1000*repmat(rho_std, 1, 1, Ndraw).*randn(size(radar.age));
accum_dt = 0.02*(1000*repmat(rho_mod, 1, 1, Ndraw) + noise_rho);

% accum_yr = (yr_top:-1:yr_end)';
% accum = zeros(length(accum_yr), size(radar.age, 2), Ndraw);
% for i = 1:size(accum, 2)
%     for j = 1:Ndraw
%         years_j = radar.age(:,i,j);
%         yr_idx = logical([diff(floor(years_j)); 0]);
%         yr_loc = find(yr_idx);
%         if length(yr_loc) - 1 > length(accum_yr)
%             n_length = length(accum_yr);
%         else
%             n_length = length(yr_loc) - 1;
%         end
%         %         yr_all = round(years_i(yr_idx(:,j),j));
%         %         yr_num = yr_all(2:end);
%         % accum_yr(1:length(yr_num),i,j) = yr_num;
%         %         accum_j = zeros(length(yr_num), 1);
%         accum_j = zeros(size(accum_yr));
%         for n = 1:n_length
%             accum_j(n) = sum(accum_dt(yr_loc(n)+1:yr_loc(n+1),i,j));
%         end
%         %         accum(1:length(yr_num),i,j) = accum_j;
%         accum(:,i,j) = accum_j(1:length(accum_yr));
%     end
% end

% Define age-depth profile as the mean of all MC age profiles for each
% trace (avoids integer year jumps in accumulation estimates and
% non-monotonically decreasing ages)
% ages = repmat(median(radar.age, 3), 1, 1, Ndraw);
ages = radar.age;

% accum_yr_init = (yr_top:-1:yr_end)';
% accum = zeros(length(accum_yr_init), size(radar.age, 2), Ndraw);
% for i = 1:size(accum, 2)
%     for j = 1:Ndraw
%         years_j = ages(:,i,j);
%         yr_idx = logical([diff(floor(years_j)); 0]);
%         yr_loc = find(yr_idx);
%         if length(yr_loc) - 1 > length(accum_yr_init)
%             n_length = length(accum_yr_init);
%         else
%             n_length = length(yr_loc) - 1;
%         end
%         %         yr_all = round(years_i(yr_idx(:,j),j));
%         %         yr_num = yr_all(2:end);
%         % accum_yr(1:length(yr_num),i,j) = yr_num;
%         %         accum_j = zeros(length(yr_num), 1);
%         accum_j = zeros(size(accum_yr_init));
%         for n = 1:n_length
%             accum_j(n) = sum(accum_dt(yr_loc(n)+1:yr_loc(n+1),i,j));
%         end
%         %         accum(1:length(yr_num),i,j) = accum_j;
%         accum(:,i,j) = accum_j(1:length(accum_yr_init));
%     end
% end

% Define initial accumulation year vector (will be iteratively modified at
% each trace in for loop)
accum_yr_init = (yr_top-1:-1:yr_end)';

% Preallocation of cell arrays for accumulation years and annual
% accumulation rate
accum_yr = cell(1, size(radar.age, 2));
accum = cell(1, size(radar.age, 2));
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
            round(1*(mean(diff(yr_loc))/12)*randn(length(yr_loc)-2, 1));
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

% %% Calculate annual accumulation for the spatially weighted composite core
% 
% % Calculate SWE accumulation at each depth interval in the weighted
% % composite core
% core_accum_dt = 0.02*(1000*core.rho);
% 
% % Find indices of integer ages within core age profile
% yr_top = floor(core.age(2));
% yr_end = ceil(core.age(end));
% core_yr = (yr_top:-1:yr_end)';
% core_yr_idx = logical([1; diff(floor(core.age))]);
% yr_loc = find(core_yr_idx);
% 
% core_accum = zeros(length(core_yr), Ndraw);
% for j = 1:Ndraw
%     
%     % Add noise to integer age locations due to uncertainty in exact point 
%     % in time of the accumulation peak, using a std dev of 1 month
%     yr_loc_j = yr_loc;
%     yr_loc_j(2:end-1) = yr_loc(2:end-1) + ...
%         round(1*(mean(diff(yr_loc))/12)*randn(length(yr_loc)-2, 1));
%     loc_idx = yr_loc_j<1;
%     yr_loc_j(loc_idx) = yr_loc(loc_idx);
%     
%     % Integrate accumulation at each depth point for each whole year in
%     % firn core
%     core_accum_j = zeros(length(core_yr), 1);
%     for n = 1:length(core_yr)
%         core_accum_j(n) = sum(core_accum_dt(yr_loc_j(n)+1:yr_loc_j(n+1)));
%     end
%     
%     % Output accumulatio results to preallocated array
%     core_accum(:,j) = core_accum_j;
% end
% 
% 
% 
% % % Calculate annual accumulation for the spatially weighted composite core
% % core_accum_dt = 0.02*(1000*core.rho);
% % core_yr_idx = logical([1; diff(floor(core.age))]);
% % yr_loc = find(core_yr_idx);
% % yr_all = round(core.age(yr_loc));
% % core_yr = yr_all(2:end);
% % core_accum = zeros(length(core_yr), 1);
% % for n = 1:length(core_yr)
% %     core_accum(n) = sum(core_accum_dt(yr_loc(n)+1:yr_loc(n+1)));
% % end

%% Output results to radar and core structures

radar.SMB_yr = accum_yr;
radar.SMB = accum;
% core.SMB_yr = core_yr;
% core.SMB = core_accum;

% SMB = struct('Easting', Easting_stack, 'Northing', Northing_stack, ...
%     'radar_accum', accum_stack, 'radar_yr', accum_yr, 'radar_ERR', ...
%     ERR_stack, 'core_accum', core_accum, 'core_yr', core_yr);
end