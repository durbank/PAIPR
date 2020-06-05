% Simplified script to import, format, and concatenate OIB .nc echograms
% into 30 km subdomains

function [files, start_idx, end_idx] = OIB_chunk(echo_dir)

% Get list of .nc files within echogram directory
wild = '*.nc';
files = dir(fullfile(echo_dir, wild));


dist_raw = 5000; %Length of single raw echogram (meters)
dist_final = 30000; %Desired length of final processing subdomain (meters)


%%% Needs some lines to determine where breaks in flight continuity occur,
%%% prior to chunking into 30 km sections


start_idx = 1:(dist_final/dist_raw):length(files);
end_idx = (dist_final/dist_raw):(dist_final/dist_raw):length(files);

%%% Later can add additional considerstions for what to do with leftover
%%% raw echogram files (total length not divisible by 6)
    