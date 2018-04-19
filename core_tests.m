% Script to generate radar accumulation data for test data for each SEAT
% firn core site

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/Data/';
        addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Addons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

% Add OIB scripts to path
addpath cresis-L1B-matlab-readers/

% Number of Monte Carlo simulations
Ndraw = 100;

% Import firn core data
[cores] = import_cores(strcat(data_path, 'Ice-cores/SEAT_cores/', ...
    'DGK_core_data.xlsx'), Ndraw);

%%

site_names = {'SEAT10_1' 'SEAT10_3' 'SEAT10_4' 'SEAT10_5' 'SEAt10_6' ...
    'SEAT11_1' 'SEAT11_2' 'SEAT11_3' 'SEAT11_4' 'SEAT11_6' 'SEAT11_7' 'SEAT11_8'};
radar_dir = strcat(data_path, 'radar/SEAT_Traverses/core-site_tests/');

for i = site_names


data_name = char(i);
file_name = strcat('layers_ku_band_',data_name, '.mat');
file = strcat(radar_dir, file_name);

% Calculate radar ages and associated other data
[radar] = radar_age(file, cores, Ndraw);

% Calculate annual accumulation rates from data
[radar] = calc_SWE(radar, Ndraw);

% Remove first/last 10 traces in radar (addresses some edge effect problems
% present in many data sets
edge = 25;

radar.Easting = radar.Easting(edge:end-edge);
radar.Northing = radar.Northing(edge:end-edge);
radar.dist = radar.dist(edge:end-edge);
radar.data_stack = radar.data_stack(:,edge:end-edge);
radar.rho_coeff = radar.rho_coeff(:,edge:end-edge);
radar.rho_var = radar.rho_var(:,edge:end-edge);
radar.data_smooth = radar.data_smooth(:,edge:end-edge);
radar.layer_vals = radar.layer_vals(:,edge:end-edge);
radar.likelihood = radar.likelihood(:,edge:end-edge);
radar.ages = radar.ages(:,edge:end-edge,:);
radar.SMB_yr = radar.SMB_yr(edge:end-edge);
radar.SMB = radar.SMB(edge:end-edge);

output_path = strcat(radar_dir, data_name,  '_radar.mat');
save(output_path, '-struct', 'radar', '-v7.3')

clearvars -except data_path Ndraw cores site_names radar_dir i

end

