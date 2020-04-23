% A function to save the output of radar processing within a parellel for
% loop (so as to preserve variable transparency). This saves both the .mat
% results of PAIPR as well as the .csv results of gamma-fitting

function [success] = parsave_all(mdata, mat_output, csv_output)


%% Fit gamma distributions to results and save to disk

% Convert Northing, Easting to lat/lon
[lat, lon] = ps2ll(mdata.Easting, mdata.Northing);

% Fit gamma distribution parameters to each year of accumulation data for
% each trace in data
table_cell = cellfun(@(lat, lon, elev, time, SMB, year) ...
    accum_distGamma(lat, lon, elev, time, SMB, year), ...
    num2cell(lat), num2cell(lon), num2cell(mdata.elev), ...
    num2cell(mdata.collect_time), mdata.SMB, mdata.SMB_yr, ...
    'UniformOutput', false);

% Concatenate results into a single table in long format
long_table = vertcat(table_cell{:});

% Create table for QC values of echogram image
QC_T = table(mdata.QC_flag, mdata.QC_med, mdata.QC_val, mdata.QC_depth_idx,...
    'VariableNames', {'QC_flag', 'QC_med', 'QC_val', 'QC_depth'});

% Add variables for QC values of echogram image
full_table = [long_table repmat(QC_T, height(long_table), 1)];

% Save the file
writetable(full_table, csv_output);

%%

save(mat_output, '-struct', 'mdata', '-v7.3')
success = true;

end