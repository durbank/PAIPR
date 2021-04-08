% Function to calculate expected value and variance of annual SMB
% distributions using a Gaussian mixture model with variable number of
% nodes (determined via peak finding)
% The function also aggregates results (by bin size) into more broad
% distributions

function [accum_table] = accum_GMix(radar, bin_size)

% Determine number of elements in each bin based on distance spacing and
% requested bin_size (bin_size in meters)
dist_interval = mean(diff(radar.dist));
stack_step = round(bin_size / dist_interval);

% Determine starting and stopping indices for different aggregated bins
start_idx = 1:stack_step:length(radar.SMB);
end_idx = [start_idx(2:end)-1 length(radar.SMB)];

% Perform gaussian mixture model fitting for each binned location
accum_table = cell(1,length(start_idx));
for i=1:length(start_idx)
    
    % Extract current bin data
    SMB_i = radar.SMB(start_idx(i):end_idx(i));
    yr_i = radar.SMB_yr(start_idx(i):end_idx(i));
    
    % Convert SMB cells to matrix
    cell_len = cellfun(@length, yr_i);
    [~,max_idx] = max(cell_len);
    SMB_mat = nan(max(cell_len), length(SMB_i),size(SMB_i{1},2));
    for j=1:length(SMB_i)
        SMB_mat(1:cell_len(j),j,:) = SMB_i{j};
    end
    
    % Stack results from different MC simulations
    SMB_mat = reshape(...
        SMB_mat, size(SMB_mat,1), size(SMB_mat,2)*size(SMB_mat,3));
    
    % Keep only results where at least 75% of simulations have coverage
    nan_sum = sum(isnan(SMB_mat),2) / size(SMB_mat,2);
    if any(nan_sum > 0.25, 'all')
        cut_idx = find(nan_sum > 0.25, 1, 'first') - 1;
        SMB_mat = SMB_mat(1:cut_idx,:);
    else
        cut_idx = size(SMB_mat, 1);
    end
    
    % Preallocation
    SMB_mu = cell(size(SMB_mat,1),1);
    SMB_std = cell(size(SMB_mat,1),1);
    SMB_scale = cell(size(SMB_mat,1),1);
    
    for j=1:length(SMB_mu)
        
        try
            
            % Remove NaNs from current data
            data_j = SMB_mat(j,:);
            data_j = data_j(~isnan(data_j));
            
            % Estimate distribution shape (cutoff at 99 percentile)
            [f,~] = ksdensity(data_j, ...
                0:round(quantile(data_j, 0.99)));
            
            % Find all peaks in distribution at least 15% the height of the
            % max peak
            p_array = findpeaks(f, 'MinPeakHeight', 0.15*max(f));
            % Another option here is to set all ksdensity values less than
            % 0.15*max(f) to zero, thus eliminating the long tails...
            
            
%             if length(p_array) < 2
%                 p_array = ones(1,2);
%             end
%
%             % Fit Gaussian mixture model using the determined number of
%             % nodes
%             GMmod = fitgmdist(data_j', length(p_array), ...
%                 'CovarianceType', 'diagonal');
%
%             %%%% Diagnostics
%             data_cells = cell(1, length(p_array));
%             for k=1:length(data_cells)
%                 data_cells{k} = normrnd(GMmod.mu(k), ...
%                     sqrt(GMmod.Sigma(k)), 1, ...
%                     round(GMmod.ComponentProportion(k)*size(SMB_mat,2)));
%             end
%             data_plt = horzcat(data_cells{:});
%
%             figure
%             hold on
%             ksdensity(data_j)
%             ksdensity(data_plt)
%             hold off
%             %%%% Diagnostics
%
%             % Take only components representing at least 10% of results
%             keep_idx = GMmod.ComponentProportion >= 0.10;
%             SMB_mu{j} = GMmod.mu(keep_idx);
%             SMB_std{j} = sqrt(GMmod.Sigma(keep_idx));
%             SMB_scale{j} = squeeze(GMmod.ComponentProportion(keep_idx));
            
            % Deterimine if mixture model is needed
            if length(p_array) > 1
                
                % Fit Gaussian mixture model using the determined number of
                % nodes
                GMmod = fitgmdist(data_j', length(p_array), ...
                    'CovarianceType', 'diagonal');
                
%                 %%%% Diagnostics
%                 data_cells = cell(1, length(p_array));
%                 for k=1:length(data_cells)
%                     data_cells{k} = normrnd(GMmod.mu(k), ...
%                         sqrt(GMmod.Sigma(k)), 1, ...
%                         round(GMmod.ComponentProportion(k) * ...
%                         size(SMB_mat,2)));
%                 end
%                 data_plt = horzcat(data_cells{:});
%                 
%                 figure
%                 hold on
%                 ksdensity(data_j)
%                 ksdensity(data_plt)
%                 hold off
%                 %%%% Diagnostics
                
                % Take only components representing at least 10% of results
                keep_idx = GMmod.ComponentProportion >= 0.10;
                SMB_mu{j} = GMmod.mu(keep_idx);
                SMB_std{j} = sqrt(GMmod.Sigma(keep_idx));
                SMB_scale{j} = squeeze(...
                    GMmod.ComponentProportion(keep_idx));
                
            else
                
                % If only one peak present, fall back to standard uni-modal
                % Gaussian model
                SMB_mu{j} = nanmean(SMB_mat(j,:));
                SMB_std{j} = nanstd(SMB_mat(j,:));
                SMB_scale{j} = 1;
                
            end
            
        catch
            
            % If mixture model fails, default to simple Gaussian
            % distribution statistics
            SMB_mu{j} = nanmean(SMB_mat(j,:));
            SMB_std{j} = nanstd(SMB_mat(j,:));
            SMB_scale{j} = 1;
            
        end
        
    end
    
    % For now, keep only dominant component (other info could be useful in
    % the future...
    [~, SMB_idx] = cellfun(@max, SMB_scale);
    mu1 = zeros(length(SMB_idx), 1);
    std1 = zeros(length(SMB_idx), 1);
    for j=1:length(SMB_idx)
        mu1(j) = SMB_mu{j}(SMB_idx(j));
        std1(j) = SMB_std{j}(SMB_idx(j));
    end
    
    % Determine average bin values for remaining variables of interest
    yr_len = size(SMB_mat,1);
    time = mean(radar.collect_time(start_idx(i):end_idx(i)));
    lat = mean(radar.Lat(start_idx(i):end_idx(i)));
    lon = mean(radar.Lon(start_idx(i):end_idx(i)));
    elev = mean(radar.elev(start_idx(i):end_idx(i)));
    
    % Assign current bin results to output table
    accum_table{i} = table(repmat(time, yr_len,1), repmat(lat, yr_len,1), ...
        repmat(lon, yr_len,1), repmat(elev, yr_len,1), ...
        yr_i{max_idx}(1:cut_idx), mu1, std1, ...
        'VariableNames', ...
        {'collect_time', 'Lat', 'Lon', 'elev', ...
        'Year', 'accum_mu', 'accum_std'});
    
end
