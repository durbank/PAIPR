% Script to import and format Operation IceBridge radar data for use with
% accum-radar scripts

addpath cresis-L1B-matlab-readers/

file = '/Volumes/WARP/Research/Antarctica/Data/IceBridge/Snow Radar/2011/IRSNO1B_20111109_02_180.nc';

mdata = load_L1B(file);