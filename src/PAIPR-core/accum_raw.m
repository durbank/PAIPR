% Blah blah blah
function [accum_table] = accum_raw(radar)

% Preallocate cell for results
accum_table = cell(1,length(radar.SMB));

for i=1:length(radar.SMB)
    
    tbl_tmp = cell2table(...
        {radar.collect_time(i), radar.Lon(i), radar.Lat(i), ...
        radar.elev(i), radar.SMB_yr{i}, num2cell(radar.SMB{i},2)}, ...
        'VariableNames', ...
        {'collect_time', 'Lon', 'Lat', 'elev', 'Year', 'accum'});
    
    j_max = cellfun(@length, tbl_tmp.Year);
    k_max = max(cellfun(@(x) size(x,2), tbl_tmp.accum{:}));
    accum_i = zeros(j_max*k_max,1);
    year_i = zeros(j_max*k_max,1);
    
    for j=1:j_max
        
        year_i((j-1)*k_max+1:(j-1)*k_max+k_max) = repmat(...
            tbl_tmp.Year{1}(j), k_max,1);
    	accum_i((j-1)*k_max+1:(j-1)*k_max+k_max) = (tbl_tmp.accum{1}{j})';
        
%         for k=1:k_max
%             
%             accum_i((j-1)*k_max+1:(j-1)*k_max+k_max) = repmat(...
%                 tbl_tmp.Year{1}(j), k_max,1);

    end
    
    accum_table{i} = table(...
        repmat(tbl_tmp.collect_time, length(accum_i),1), ...
        cast(repmat(tbl_tmp.Lon, length(accum_i),1), 'single'), ...
        cast(repmat(tbl_tmp.Lat, length(accum_i),1), 'single'), ...
        cast(repmat(tbl_tmp.elev, length(accum_i),1), 'single'), ...
        cast(year_i, 'int16'), cast(accum_i, 'single'), ...
        'VariableNames', tbl_tmp.Properties.VariableNames);
    
end

end