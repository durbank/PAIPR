function [accum_table] = accum_GMix(radar, bin_size)

dist_interval = mean(diff(radar.dist));
stack_step = round(bin_size / dist_interval);

start_idx = 1:stack_step:length(radar.SMB);
end_idx = [start_idx(2:end)-1 length(radar.SMB)];



% Perform gaussian mixture model fitting for each binned location
accum_table = cell(1,length(start_idx));
for i=1:length(start_idx)
    
    SMB_i = radar.SMB(start_idx(i):end_idx(i));
    yr_i = radar.SMB_yr(start_idx(i):end_idx(i));
    
    
    cell_len = cellfun(@length, yr_i);
    [~,max_idx] = max(cell_len);
    SMB_mat = nan(max(cell_len), length(SMB_i),size(SMB_i{1},2));
    for j=1:length(SMB_i)
        SMB_mat(1:cell_len(j),j,:) = SMB_i{j};
    end
    
    SMB_mat = reshape(...
        SMB_mat, size(SMB_mat,1), size(SMB_mat,2)*size(SMB_mat,3));
    nan_sum = sum(isnan(SMB_mat),2) / size(SMB_mat,2);
    SMB_mat = SMB_mat(nan_sum <= 0.33,:);
    
    SMB_mu = cell(size(SMB_mat,1),1);
    SMB_std = cell(size(SMB_mat,1),1);
    SMB_scale = cell(size(SMB_mat,1),1);
    
    for j=1:length(SMB_mu)
        
        try
            
            [f,~] = ksdensity(SMB_mat(j,:));
            p_array = findpeaks(f, 'MinPeakHeight', 0.15*max(f));
            % Another option here is to set all ksdensity values less than
            % 0.15*max(f) to zero, thus eliminating the long tails...
            if length(p_array) < 2
                p_array = ones(1,2);
            end
        
        
            GMmod = fitgmdist(SMB_mat(j,:)', length(p_array), ...
                'CovarianceType', 'diagonal');
    
%             %%%% Diagnostics
%             data_cells = cell(1, length(p_array));
%             for k=1:length(data_cells)
%                 data_cells{k} = normrnd(GMmod.mu(k), ...
%                     sqrt(GMmod.Sigma(k)), 1, ...
%                     round(GMmod.ComponentProportion(k)*size(SMB_mat,2)));
%             end
%             data = horzcat(data_cells{:});
%         
%             figure
%             hold on
%             ksdensity(SMB_mat(j,:))
%             ksdensity(data)
%             hold off
%             %%%% Diagnostics

            % Take only components representing at least 10% of results
            keep_idx = GMmod.ComponentProportion >= 0.10;
            SMB_mu{j} = GMmod.mu(keep_idx);
            SMB_std{j} = sqrt(GMmod.Sigma(keep_idx));
            SMB_scale{j} = squeeze(GMmod.ComponentProportion(keep_idx));
            
        catch
            
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
    
    yr_len = size(SMB_mat,1);
    time = mean(radar.collect_time(start_idx(i):end_idx(i)));
    lat = mean(radar.Lat(start_idx(i):end_idx(i)));
    lon = mean(radar.Lon(start_idx(i):end_idx(i)));
    elev = mean(radar.elev(start_idx(i):end_idx(i)));
    
    
    accum_table{i} = table(repmat(time, yr_len,1), repmat(lat, yr_len,1), ...
        repmat(lon, yr_len,1), repmat(elev, yr_len,1), ...
        yr_i{max_idx}(nan_sum <= 0.33), mu1, std1, ...
        'VariableNames', ...
        {'collect_time', 'Lat', 'Lon', 'elev', ...
        'Year', 'accum_mu', 'accum_std'}); 
    
end
