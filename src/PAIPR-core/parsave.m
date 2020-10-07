% A function to save the output of radar processing within a parellel for
% loop (so as to preserve variable transparency). This saves both the .mat
% results of PAIPR as well as the .csv results of gamma-fitting

function [success] = parsave(mdata, csv_output, varargin)

%% Fit gamma distributions to results and save to disk

distribution = 'mixture';

switch distribution

    case 'gamma'
        % Fit gamma distribution parameters to each year of accumulation 
        % data for each trace in data
        table_cell = cellfun(@(lat, lon, elev, time, SMB, year) ...
            accum_distGamma(lat, lon, elev, time, SMB, year), ...
            num2cell(lat), num2cell(lon), num2cell(mdata.elev), ...
            num2cell(mdata.collect_time), mdata.SMB, mdata.SMB_yr, ...
            'UniformOutput', false);

        % Concatenate results into a single table in long format
        long_table = vertcat(table_cell{:});
        
    case 'gaussian'
        % Fit gaussian distributions to each year of accumulation data
        % (also aggregates traces by given distance bin in meters)
        bin_size = 200;
        table_cell = accum_Gauss(mdata, bin_size);
        
        % Converts data to a long form space-time table
        long_table = vertcat(table_cell{:});
        
    case 'mixture'
        % Fit Gaussian mixture model distributions to each year of
        % accumulation data (also aggregates traces by given distance bins)
        bin_size = 200;
        table_cell = accum_GMix(mdata, bin_size);
        
        % Converts data to a long form space-time table
        long_table = vertcat(table_cell{:});
end

% Create table for QC values of echogram image
QC_T = table(mdata.QC_flag, mdata.QC_med, mdata.QC_val, ...
    mdata.QC_yr, 'VariableNames', ...
    {'QC_flag', 'QC_med', 'QC_val', 'QC_yr'});

% Add variables for QC values of echogram image
full_table = [long_table repmat(QC_T, height(long_table), 1)];

% Save the file
writetable(full_table, csv_output);

%%

if length(varargin) == 1
    mat_output = varargin{1};
    save(mat_output, '-struct', 'mdata', '-v7.3')
end


success = true;

end