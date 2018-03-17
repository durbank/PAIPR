function [radar, core] = calc_SWE(radar, core)

% Generate mean age-depth profiles for each trace from MC simulations
ages = mean(radar.age, 3);

% Define first year with complete accumulation data and last year with
% robust accumulation data
yr_top = floor(max(ages(1,:)));
yr_end = ceil(min(ages(end,:)));

% Generate density with depth model for radar data
rho_ice = 0.917;
rho_0 = mean(core.rho(1:10));
rho_coeff = coeffvalues(radar.rho_fit);
rho_mod = rho_coeff(1)*radar.depth.^rho_coeff(2) + rho_coeff(3);
% rho_mod = rho_ice - (rho_ice-rho_0)*exp(rho_coeff.*radar.depth_interp);

% Generate variance model for density with depth, and calculate standard
% deviation with depth for density model
rho_res = rho_coeff(1)*core.depth.^rho_coeff(2) + rho_coeff(3);
% rho_res = core.rho - (rho_ice - (rho_ice-rho_0)*exp(rho_coeff*core.depth));
rho_var = movvar(rho_res, round(length(rho_res)/3));
fun = @(a, x) (max(rho_var)-mean(rho_var(end-50:end)))./(a.^x) + ...
    mean(rho_var(end-50:end));
var_fit = fit(core.depth, rho_var, fun, 'StartPoint', 1.5);
var_coeff = coeffvalues(var_fit);
std_mod = 1000*sqrt(abs((max(rho_var)-mean(rho_var(end-50:end)))./...
    (var_coeff.^radar.depth) + mean(rho_var(end-50:end))));

% % Calculate accumulation at each depth interval (0.02 m) with simulated
% % noise
noise = repmat(std_mod, 1, size(radar.age, 2)).*randn(size(radar.age));
accum_dt = 0.02*(1000*repmat(rho_mod, 1, size(ages, 2)) + noise);
% accum_dt = 0.02*(1000*repmat(rho_mod, 1, size(radar.age, 2), Ndraw) + noise);

accum_yr = (yr_top:-1:yr_end)';
accum = zeros(length(accum_yr), size(radar.age, 2));
for i = 1:size(accum, 2)
    years_i = ages(:,i);
    yr_idx = logical([diff(floor(years_i)); 0]);
    yr_loc = find(yr_idx);
    if length(yr_loc) - 1 > length(accum_yr)
        n_length = length(accum_yr);
    else
        n_length = length(yr_loc) - 1;
    end
    %         yr_all = round(years_i(yr_idx(:,j),j));
    %         yr_num = yr_all(2:end);
    % accum_yr(1:length(yr_num),i,j) = yr_num;
    %         accum_j = zeros(length(yr_num), 1);
    accum_i = zeros(size(accum_yr));
    for n = 1:n_length
        accum_i(n) = sum(accum_dt(yr_loc(n)+1:yr_loc(n+1),i));
    end
    %         accum(1:length(yr_num),i,j) = accum_j;
    accum(:,i) = accum_i(1:length(accum_yr));
    
end

% Truncate accum to years where all traces have data for all years
accum_yr = accum_yr(all(accum, 2));
if max(diff(accum_yr)) > 1
    disp('Issue with accumulation year trunctation')
end
accum = accum(all(accum ,2),:);

% % Stack adjacent records based on the lateral averaging window size
% window = 50;
% stack_idx = 1:window:size(accum, 2);
% accum_stack = zeros(length(accum_yr), length(stack_idx)-1);
% ERR_stack = zeros(length(accum_yr), length(stack_idx)-1);
% Northing_stack = zeros(1, length(stack_idx)-1);
% Easting_stack = zeros(1, length(stack_idx)-1);
% for i = 1:length(stack_idx)-1
%     accum_stack(:,i) = mean(accum(:,stack_idx(i):stack_idx(i+1)), 2);
%     ERR_stack(:,i) = 1.96*std(accum(:,stack_idx(i):stack_idx(i+1)), [], 2)/sqrt(window);
%     Northing_stack(i) = mean(radar.Northing(stack_idx(i):stack_idx(i+1)));
%     Easting_stack(i) = mean(radar.Easting(stack_idx(i):stack_idx(i+1)));
% end

% Calculate annual accumulation for the spatially weighted composite core
core_accum_dt = 0.02*(1000*core.rho);
core_yr_idx = logical([1; diff(floor(core.age))]);
yr_loc = find(core_yr_idx);
yr_all = round(core.age(yr_loc));
core_yr = yr_all(2:end);
core_accum = zeros(length(core_yr), 1);
for n = 1:length(core_yr)
    core_accum(n) = sum(core_accum_dt(yr_loc(n)+1:yr_loc(n+1)));
end

radar.SMB_yr = accum_yr;
radar.SMB = accum;
core.SMB_yr = core_yr;
core.SMB = core_accum;

% SMB = struct('Easting', Easting_stack, 'Northing', Northing_stack, ...
%     'radar_accum', accum_stack, 'radar_yr', accum_yr, 'radar_ERR', ...
%     ERR_stack, 'core_accum', core_accum, 'core_yr', core_yr);
end