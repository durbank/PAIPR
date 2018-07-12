% Script used to generate processed age-depth and accumulation estimates
% from OIB 2011 snow radar between SEAT10-4 and SEAT10-6

% Load OIB data from file (data used was previously generated using
% OIB_concat.m)
OIB_file = '/media/durbank/WARP/Research/Antarctica/Data/IceBridge/SEAT10_4to10_6/2011_SNO_all.mat';
radar = load(OIB_file);

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = '/media/durbank/WARP/Research/Antarctica/Data/Ice-cores/SEAT_cores/SEAT_cores.mat';
cores = load(core_file);

% Define number of MC simulations
Ndraw = 100;

% Generate estimates of annual horizons and age-depth scale estimates
[radar] = OIB_age(radar, cores, Ndraw);

% Calculate annual accumulation rates from data
[radar] = calc_SWE(radar, Ndraw);
