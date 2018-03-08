% Script for testing the import and processing of OIB data with my codes

PC_true = ispc;
switch PC_true
    case true
        data_path = 'D:/Research/Antarctica/WAIS Variability/';
        addon_path = 'D:/Research/Antarctica/WAIS Variability/Addons/';
    case false
        data_path = '/Volumes/WARP/Research/Antarctica/WAIS Variability/';
        addon_path = '/Users/Durbank/Documents/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

addpath cresis-L1B-matlab-readers/

% Files IRSNO1B_20111109_02_211 through ...02_227 roughly follows the
% SEAT10-1 to SEAT10-6 transect. The same holds true for the following:
%   IRSNO1B_20161109_02_320 through IRSNO1B_20161109_02_336
%   Ku-band: IRKUB1B_20161109_02_320 through IRKUB1B_20161109_02_336

% Files IRSNO1B_20111109_02_242 through ...02_272 exactly follows the 
% SEAT10-4 - SEAT10-6 transect. The same holds true for the following:
%   IRSNO1B_20161109_02_350 through IRSNO1B_20161109_02_381
%   Ku-band: IRKUB1B_20161109_02_350 through IRKUB1B_20161109_02_381

% Path of the OIB file to process
file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow radar/2011/IRSNO1B_20111109_02_211.nc';

% Number of simulations to perform on age-depth Monte Carlo
Ndraw = 100;

%%

% Import firn core data
[cores] = import_cores(strcat(data_path, ['SEAT_cores' filesep ...
    'DGK_core_data.xlsx']));

% Import and process OIB data
[radar] = OIB_import(file);




