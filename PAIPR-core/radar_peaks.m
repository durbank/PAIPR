

function [peak_prom, peak_width] = radar_peaks(radar, vert_res)

%% Find depth, width, and prominence of peaks for each radar trace

% Preallocate arrays for various components
peak_prom = zeros(size(radar.data_smooth));
peak_width = zeros(size(radar.data_smooth));
% Proms = cell(1, size(radar.data_smooth, 2));
% widths = cell(1,size(radar.data_smooth, 2));
% depths = cell(1, size(radar.data_smooth, 2));
% depth_idx = cell(1, size(radar.data_smooth, 2));

for i = 1:size(radar.data_smooth, 2)
    data_i = radar.data_smooth(:,i);
    
    % Prominence threshold for peaks
    minProm = 0.25;
    
    % Min distance between peaks (in meters)
    minDist = 0.08;
    
    % Find peak statistics in each trace based on criteria
    [~, peaks_idx_i, widths_i, Prom_i] = findpeaks(data_i, ...
        'MinPeakProminence', minProm, ...
        'MinPeakDistance', minDist/vert_res, 'WidthReference', 'halfprom');

    % Add peak prominence and width values to relevent matrices
    peak_prom(peaks_idx_i,i) = Prom_i;
    peak_width(peaks_idx_i,i) = widths_i;
    
%     % Add values to relevent cell arrays
%     Proms{i} = Prom_i;
%     widths{i} = widths_i;
%     depths{i} = radar.depth(peaks_idx_i);
%     depth_idx{i} = peaks_idx_i;
end


end