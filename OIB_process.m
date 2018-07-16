% Script used to generate processed age-depth and accumulation estimates
% from OIB 2011 snow radar between SEAT10-4 and SEAT10-6

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'work';
%         computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
                
            case 'laptop'
                data_path = 'C:/Users/durba/Documents/Research/Antarctica/Data/';
                addon_path = 'C:/Users/durba/Documents/MATLAB/Addons/';
        end
        
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = fullfile(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))
% Add export_fig to path
addon_folder = fullfile(addon_path, 'altmany-export_fig-cafc7c5/');
addpath(genpath(addon_folder))

% Load OIB data from file (data used was previously generated using
% OIB_concat.m)
OIB_file = fullfile(data_path, 'IceBridge/SEAT10_4to10_6/2011_SNO_all.mat');
radar = load(OIB_file);

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Define number of MC simulations
Ndraw = 100;

% Generate estimates of annual horizons and age-depth scale estimates
[radar] = OIB_age(radar, cores, Ndraw);

% Calculate annual accumulation rates from data
[radar] = calc_SWE(radar, Ndraw);
