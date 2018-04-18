% Function to import core data from .xlsx file (specific format required)
% and arrange it in a structure for future MATLAB analysis

% core_file = '/Volumes/WARP DRIVE/Research/Antarctica/WAIS Variability/SEAT_cores/core_data_interp.xlsx';
function [cores] = import_cores(core_file, Ndraw)

% % Add Antarctic Mapping Toolbox (AMT) to path
% addon_path = ['Addons' filesep 'AntarcticMappingTools_v5.03' filesep];
% addpath(genpath(addon_path))

% Import basic data about firn cores in .xlsx file (name, location,
% elevation, etc) and places inside a data structure
[data, name] = xlsread(core_file, 'Site Information');
name = name(1,2:end);
cores = struct('name', {name}, 'lat', data(1,:), 'lon', data(2,:), 'elev', ...
    data(3,:));

% Convert core lat/lon locations to polar sterographic projection (meters)
[cores.Easting, cores.Northing] = ll2ps(cores.lat, cores.lon);

% Import additional core information (depth, age, density, isotopes, etc), 
% fill in NaN/missing values in density and isotope data, and place data in 
% the 'cores' data structure
for i = 1:numel(cores.name)
    [data, ~] = xlsread(core_file, cores.name{i});
    cores.(cores.name{i}) = struct('name', cores.name(i), 'Easting', ...
        cores.Easting(i), 'Northing', cores.Northing(i), 'depth', ...
        data(:,1), 'age', data(:,2), 'rho', fillmissing(data(:,3), 'pchip'),...
        'dD', fillmissing(data(:,4), 'pchip'), ...
        'd18O', fillmissing(data(:,5), 'pchip'));
    
    %% Calculate annual accumulation in current core (mm per year)
    
    % Define ith core
    core_i = cores.(cores.name{i});
    
    % Calculate SWE accumulation at each depth interval in the core,
    % assuming 5% error (95% confidence interval) in density measurement
    error = 0.05;
    rho_i = repmat(core_i.rho, 1, Ndraw);
    rho_std = (error/2)*core_i.rho;
    rho_noise = repmat(rho_std, 1, Ndraw).*randn(size(rho_i));
    core_accum_dt = 0.02*(1000*(core_i.rho + rho_noise));
    
    % Add Gaussian noise to core age, assuming a std of +/- 0.5 year per
    % decade
%     duration = max(core_i.age) - min(core_i.age);
%     std_bott = (0.5/10)*duration;
%     age_std = (0:std_bott/(length(core_i.age)-1):std_bott)';
% 
%     std_time = age_std.*randn(length(core_i.age), Ndraw) + ...
%         repmat((0:duration/(length(core_i.age)-1):duration)', 1, Ndraw);
%     
%     duration = max(core_i.age) - min(core_i.age);
%     std_bott = (0.5/10)*duration;
%     age_std = (0:std_bott/(length(core_i.age)-1):std_bott)';
    

   
    ages = repmat(core_i.age, 1, Ndraw);
    
    
    % Find indices of integer ages within core age profile
    yr_top = floor(core_i.age(1));
    yr_end = ceil(min(ages(end,:)));
    core_yr_init = (yr_top-1:-1:yr_end)';
    
%     core_yr_idx = logical([1; diff(floor(core_i.age))]);
%     yr_loc = find(core_yr_idx);
    
    core_accum = zeros(length(core_yr_init), Ndraw);
    for j = 1:Ndraw
        
        years_j = ages(:,j);
        yr_idx = logical([diff(floor(years_j)); 0]);
        
        % Add noise to integer age locations due to uncertainty in exact point
        % in time of the accumulation peak, using a std dev of 1 month
        
        
        
        
        % Calculate indices of integer ages for jth simulation of the ith
        % core (with added noise from uncertainty in exact point in time
        % of the accumulation peak, using a std dev of 1 month)
        yr_loc = find(yr_idx);
        loc_temp = yr_loc;
        loc_temp(2:end-1) = yr_loc(2:end-1) + ...
            round(1*(mean(diff(yr_loc))/12)*randn(length(yr_loc)-2, 1));
        loc_idx = loc_temp<1;
        loc_temp(loc_idx) = yr_loc(loc_idx);
        yr_loc = loc_temp;
        
        % Integrate accumulation at each depth point for each whole year in
        % firn core realization
        nlength = length(yr_loc) - 1;
        core_accum_j = zeros(nlength, 1);
        for n = 1:nlength
            core_accum_j(n) = sum(core_accum_dt(yr_loc(n)+1:yr_loc(n+1),j));
        end
        
        % Output accumulatio results to preallocated array
        core_accum(1:nlength,j) = core_accum_j;
    end
    
    % Output results to respective arrays
    accum_idx = find(all(core_accum, 2), 1, 'last');
    accum_clip = core_accum(1:accum_idx,:);
    core_yr = core_yr_init(1:accum_idx);
    core_accum = accum_clip;
    
    cores.(cores.name{i}).SMB_yr = core_yr;
    cores.(cores.name{i}).SMB = core_accum;
    
end

end