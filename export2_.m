% Function to format and prepare for export PAIPR-processed data to various
% different forms

function [data_struct] = export2_(input_dir, varargin)

% varargin = {'Easting', 'Northing', 'SMB', 'SMB_yr'};

% % GUI to specify input data directory
% input_dir = uigetdir('../', ...
%     'Select directory from which to input PAIPR-processed data');

% Find file names for previously processed PAIPR data
wild = '*.mat';
files = dir(fullfile(input_dir, wild));

% Preallocate structure to house data
data_struct = struct();

for i = 1:length(varargin)
    
    arg_name = varargin{i};
    
    try
        tmp_load = load(fullfile(files(1).folder, files(1).name), arg_name);
        tmp_data = tmp_load.(arg_name);
        
        switch class(tmp_data)
            
            case 'double'
                % Preallocate array of sufficient size for data
                pointer_vect = false(1, length(files)*2*(50*1000/25));
                arg_data = zeros(1, length(files)*2*(50*1000/25));
                
                for j = 1:length(files)
                    
                    % Load relevent data from current data file
                    tmp_load = load(fullfile(files(j).folder, files(j).name), arg_name);
                    tmp_data = tmp_load.(arg_name);
                    
                    % Find position of last data entered into preallocated array
                    next_idx = sum(logical(pointer_vect)) + 1;
%                     next_idx = sum(find(var_data)) + 1;
                    
                    % Fill current iteration data into preallocated array
                    arg_data(next_idx:next_idx+length(tmp_data)-1) = tmp_data;
                    pointer_vect(next_idx:next_idx+length(tmp_data)-1) = tmp_data;
                end
                
                % Find and remove empty data indices (usually from excess preallocation)
                arg_data = arg_data(pointer_vect);
%                 keep_idx = find(logical(pointer_vect));
%                 var_data = var_data(keep_idx);
%                 last_idx = find(var_data, 1, 'last');
%                 var_data = var_data(1:last_idx);
                
            case 'cell'
                % Preallocate array of sufficient size for data
                arg_data = cell(1, length(files)*2*(50*1000/25));
                
                for j = 1:length(files)
                    
                    % Load relevent data from current data file
                    tmp_load = load(fullfile(files(j).folder, files(j).name), arg_name);
                    tmp_data = tmp_load.(arg_name);
                    
                    % Find position of last data entered into preallocated arrays
                    next_idx = sum(~cellfun(@isempty, arg_data)) + 1;
                    
                    % Fill current iteration data into preallocated array
                    arg_data(next_idx:next_idx+length(tmp_data)-1) = tmp_data;
                end
                
                % Find and remove empty SEAT indices (usually from excess preallocation)
                last_idx = sum(~cellfun(@isempty, arg_data));
                arg_data = arg_data(1:last_idx);
        end
        
        % Add ith data variable to the data structure
        data_struct.(arg_name) = arg_data;
        
    catch
        sprintf("Warning: Could not find/import '%s' data for files in the %s directory", ...
            arg_name, input_dir)
    end
end

end