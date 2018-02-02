

% Assigns the path to the data and add-on folders, based on whether running
% from Durban's Macbook Pro or a lab PC
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

% Directory to radar files of interest
radar_dir = strcat(data_path, ['SEAT_Traverses' filesep 'SEAT2010Kuband'...
    filesep 'ProcessedSEAT2010' filesep 'grid_SEAT10_4' filesep]);

% List all files matching 'wild' within radar directory
wild = 'layers*';
files = dir(strcat(radar_dir, wild));

% Import firn core data
[cores] = import_cores(strcat(data_path, ['SEAT_cores' filesep ...
    'DGK_core_data.xlsx']));

% Select individual radar file for analysis
i = randi(length(files));
file = strcat([files(i).folder filesep], files(i).name);

% Number of simulations to perform on age-depth Monte Carlo
Ndraw = 100;

% Calculate radar ages and associated other data
[radar, core] = age_MC(file, cores, Ndraw);

%Add radon scripts to path
addpath LTARE_codes/

% Generate binary image of annual layers
[layers, data_pts] = DGK_horizons(radar.data_smooth);








%%

% % Define start/end indices for the different blocks of data
% block_height = floor(size(radar.data_smooth, 1)/10);
% block_width = floor(size(radar.data_smooth, 2)/20);
% 
% % rows = buffer(1:size(radar.data_smooth, 1), window_height, round(0.25*window_height));
% % start_row = [1 rows(1,2:end)];
% % end_row = [rows(end,1:end-1) size(radar.data_smooth, 1)];
% cols = buffer(1:size(radar.data_smooth, 2), block_width, round(0.33*block_width));
% start_col = [1 cols(1,2:end)];
% end_col = [cols(end,1:end-1) size(radar.data_smooth, 2)];
% 
% start_row = 1:block_height:size(radar.data_smooth, 1)-block_height;
% end_row = [start_row(2:end) size(radar.data_smooth, 1)];
% % start_col = 1:window_width:size(radar.data_smooth, 2)-window_width;
% % end_col = [start_col(2:end) size(radar.data_smooth, 2)];
% 
% Prominence = cell(length(start_row), length(start_col));
% X_all = cell(length(start_row), length(start_col));
% Y_all = cell(length(start_row), length(start_col));
% Y_idx = cell(length(start_row), length(start_col));
% for i=1:length(start_row)
%     for j = 1:length(start_col)
%         data = radar.data_smooth(start_row(i):end_row(i),start_col(j):end_col(j));
%         
%         % % Plot the radargram
%         % figure('Position', [200 200 1500 800])
%         % imagesc(radar.dist, radar.depth_interp, radar.data_smooth, [-1.5 1.5])
%         % colorbar
%         
%         % Calculate the radon transform data
%         theta = 85:0.1:95;
%         [R_full, xp_full] = radon(data, theta);
%         
%         % Remove rows only containing zeroes
%         idx_zero = any(R_full, 2);
%         R_zero = R_full(idx_zero,:);
%         xp_zero = xp_full(idx_zero,:);
%         
%         % % Plot radon sinogram
%         % figure('Position', [200 200 1500 800])
%         % imagesc(theta, xp_zero, R_zero, [-100 1000])
%         % colormap(hot)
%         % colorbar
%         
%         % % Inverse transform the max radon data, and interpolate to original data
%         % % dimensions
%         % I1 = iradon(R_zero, theta);
%         % I = imresize(I1, size(data));
%         %
%         % figure('Position', [200 200 1500 800])
%         % imagesc(radar.dist, radar.depth_interp, I)
%         % colorbar
%         
%         % % Find the max values of the radon transform (row-wise)
%         % [~, col_idx] = max(abs(R_zero), [], 2);
%         % R_max = R_zero(col_idx);
%         % theta_max = (theta(col_idx))';
%         
%         % Calculate the median theta for the max R values, and find the location of
%         % the given theta values nearest to the median theta
%         [~, col_idx] = max(abs(R_zero), [], 2);
%         theta_max = median(theta(col_idx));
%         [~, max_idx] = min(abs(theta - theta_max));
%         
%         % Define R_max matrix, given the median theta of R_max
%         R_max = R_zero(:, max_idx);
%         
%         % Find the magnitude and location of peaks in the RT data
% %         minProm = quantile(R_max, 0.05) - min(R_max);
% %         [~, peakLoc, ~, peakProm] = findpeaks(R_max, 'minPeakProminence', minProm);
%         [~, peakLoc, ~, peakProm] = findpeaks(R_max);
%         xp_peaks = xp_zero(peakLoc);
%         
%         Prominence{i,j} = peakProm;
%         
% %         % Convert peak prominence to probability based on logistic function (may
% %         % need to adapt this, as both peak strength and lateral continuity
% %         % currently affect the estimate)
% %         B = quantile(zscore(peakProm), 0.99);
% %         P = 1./(1 + exp(-B*zscore(peakProm)));
% %         
% %         % Plot probability vs radon brightness
% %         figure
% %         plot(peakProm, P, 'bo')
% %         xlabel('Max radon brightness')
% %         ylabel('Assigned probability')
%         
%         % Draw layers based on RT calculations
%         block_sz = size(data);
%         Ymid = start_row(i) + ceil(block_sz(1)/2) - 1;
%         Xmid = start_col(j) + ceil(block_sz(2)/2) - 1;
%         m = tand(90+theta_max);
%         Yo = Ymid + xp_peaks*sind(theta_max);
%         Xo = Xmid + xp_peaks*cosd(theta_max);
%         b = mean(Yo - m.*Xo, 2);
%         x = start_col(j):end_col(j);
%         y = round(m*x + b);
%         y(y<start_row(i)) = NaN;
%         
%         X_all{i,j} = x;
%         Y_all{i,j} = y;
% %         Y_idx{i,j} = y>=start_row(i);
%         
%     end
% end
% 
% % Plot layers
% figure('Position', [200 200 1500 800])
% imagesc(radar.data_smooth, [-2 2])
% % imagesc([0 radar.dist(end)], [0 radar.depth_interp(end)], radar.data_smooth, [-2 2])
% % colorbar
% % xlabel('Distance along profile (m)')
% % ylabel('Depth (m)')
% hold on
% for i = 1:length(start_row)
%     for j = 1:length(start_col)
%         if ~isempty(Y_all{i,j})
%             plot(X_all{i,j}, Y_all{i,j}, 'r', 'LineWidth', 2)
%         end
%     end
% end
% set(gca, 'Ydir', 'reverse')
% hold off
% 
% 
% 
% %%% HOW TO COMBINE THE TWO UNCERTAINTIES %%%
% % P_peaks^(1-P_radon)
% % Could also look into conjugate priors (Bernoulli and Beta distributions)
% 
% % %%
% 
% % % % Plot radargram surface
% % % figure('Position', [200 200 1500 800])
% % % mesh(radar.dist, radar.depth_interp, radar.data_smooth)
% % % colorbar
% % % xlabel('Distance along profile (m)')
% % % ylabel('Depth (m)')
% % 
% % % Separate horizon estimates into groups of consecutive segments
% % layer_labeled = bwlabel(ann_horizon_idx);
% % % also use regionprops function
% 
% 
% % 
% % 
% % % Use radon transform values to reconstruct lines of annual layer guesses
% % % by taking each R_max value and spreading it evenly throughout a generated
% % % line. The generated line is produced using it's angle, col location, and
% % % xp value
% % n_col = size(radar.data_smooth, 2);


