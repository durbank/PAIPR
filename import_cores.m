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
    
    %% Calculate annual accumulation ith core
    
    % Define ith core
    core_i = cores.(cores.name{i});
    
    % Calculate SWE accumulation at each depth interval in the weighted
    % composite core
    core_accum_dt = 0.02*(1000*core_i.rho);
    
    % Find indices of integer ages within core age profile
    yr_top = floor(core_i.age(2));
    yr_end = ceil(core_i.age(end));
    core_yr = (yr_top-1:-1:yr_end)';
    core_yr_idx = logical([1; diff(floor(core_i.age))]);
    yr_loc = find(core_yr_idx);
    
    core_accum = zeros(length(core_yr), Ndraw);
    for j = 1:Ndraw
        
        % Add noise to integer age locations due to uncertainty in exact point
        % in time of the accumulation peak, using a std dev of 1 month
        yr_loc_j = yr_loc;
        yr_loc_j(2:end-1) = yr_loc(2:end-1) + ...
            round(1*(mean(diff(yr_loc))/12)*randn(length(yr_loc)-2, 1));
        loc_idx = yr_loc_j<1;
        yr_loc_j(loc_idx) = yr_loc(loc_idx);
        
        % Integrate accumulation at each depth point for each whole year in
        % firn core
        core_accum_j = zeros(length(core_yr), 1);
        for n = 1:length(core_yr)
            core_accum_j(n) = sum(core_accum_dt(yr_loc_j(n)+1:yr_loc_j(n+1)));
        end
        
        % Output accumulatio results to preallocated array
        core_accum(:,j) = core_accum_j;
    end
    
    
    
%     % Calculate annual accumulation for each core
%     core_accum_dt = 0.02*(1000*core_i.rho);
%     core_yr_idx = logical([1; diff(floor(core_i.age))]);
%     yr_loc = find(core_yr_idx);
%     yr_all = round(core_i.age(yr_loc));
%     core_yr = yr_all(2:end);
%     core_accum = zeros(length(core_yr), 1);
%     for n = 1:length(core_yr)
%         core_accum(n) = sum(core_accum_dt(yr_loc(n)+1:yr_loc(n+1)));
%     end
    
    cores.(cores.name{i}).SMB_yr = core_yr;
    cores.(cores.name{i}).SMB = core_accum;
    
end

end