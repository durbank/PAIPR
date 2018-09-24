% Wrapper function to import and process an arbitrary number of radargram
% files within a single directory

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'work';
        %         computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
                
            case 'laptop'
                data_path = 'C:/Users/durba/Documents/Research/Antarctica/Data/';
                addon_path = 'C:/Users/durba/Documents/MATLAB/Addons/';
        end
        
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = fullfile(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))
% Add export_fig to path
addon_folder = fullfile(addon_path, 'altmany-export_fig-cafc7c5/');
addpath(genpath(addon_folder))
% Add CReSIS OIB MATLAB reader functions to path
addon_folder = fullfile(addon_path, 'cresis-L1B-matlab-readers/');
addpath(genpath(addon_folder))

% Load core data from file (data used was previously generated using
% import_cores.m)
core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
cores = load(core_file);

% Define number of Monte Carlo simulations to perform
Ndraw = 100;

% Define path to the directory containing radar data (relative to the
% 'data' directory path)
radar_dir = fullfile(data_path, 'radar/SEAT_Traverses/SEAT2010Kuband/', ...
    'SEAT10_4toSEAT10_6');
% radar_dir = fullfile(data_path, 'IceBridge/SEAT10_4to10_6/2011_SNO');

%%

radar_ALL = radar_format(radar_dir);

for i = 1:length(radar_ALL)
    
    [radar_tmp] = radar_RT(radar_ALL(i).segment, cores, Ndraw);
    [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
    
end









%%
% % List all files matching 'wild' within radar directory
% wild = '*.mat';
% files = dir(fullfile(radar_dir, wild));
% 
% format = 'MATLAB';
% source = 'SEAT';
% 
% % If directory has no '.mat' files, additionally check for '.nc' files
% if isempty(files)
%     wild = '*.nc';
%     files = dir(fullfile(radar_dir, wild));
%     format = 'NetCDF';
%     source = 'OIB';
% end
% 
% %%
% 
% % Preallocate arrays for continuous lat/lon positions for all data in
% % directory
% lat = [];
% lon = [];
% 
% % Preallocate array for the radargram length (in data bins) for each
% % data file in directory
% file_length = zeros(1, length(files));
% 
% switch format
%     case 'MATLAB'
%         
%         % Iteratively load lat/lon positions for each data file, and add
%         % corresponding values to arrays
%         for i = 1:length(files)
%             data_i = load(fullfile(files(i).folder, files(i).name), ...
%                 'lat', 'lon');
%             file_length(i) = length(data_i.lat);
%             lat = [lat data_i.lat];
%             lon = [lon data_i.lon];
%         end
%         
%     case 'NetCDF'
%         
%         % Iteratively load lat/lon positions for each data file, and add
%         % corresponding values to arrays
%         for i = 1:length(files)
%             lat_i = ncread(fullfile(files(i).folder, files(i).name), 'lat');
%             lon_i = ncread(fullfile(files(i).folder, files(i).name), 'lon');
%             file_length(i) = length(lat_i);
%             lat = [lat lat_i'];
%             lon = [lon lon_i'];
%         end
%         
% end
% 
% % Calculate the cummulative radargram length (in data bins) for the data
% % files in directory
% file_end_idx = cumsum(file_length);
% 
% % Replace data without valid location values (outside Antarctic Circle)
% % with nearest preceding valid location (these missing data will later be
% % removed in processing)
% invld_idx = lat>=-65 | lat<-90;
% lat(invld_idx) = NaN;
% lon(invld_idx) = NaN;
% lat = fillmissing(lat, 'previous');
% lon = fillmissing(lon, 'previous');
% 
% % Calucate distance along traverse (in meters)
% d = pathdist(lat, lon);
% 
% % Find indices to break radargrams based on absence of data across an
% % extended distance (greater than 500 m)
% break_idx = [0 find(diff(d)>500) length(lat)];
% 
% % Set the minimum length needed for radargram processing and radargram
% % overlap interval (in meters)
% length_min = 30000;
% overlap = 5000;
% 
% % Preallocate and initialize variables for breakpoint calculation
% breaks_new = cell(1, length(break_idx)-1);
% j = 1;
% 
% % For each continuous section of radargram (no significant breaks),
% % determine additional breakpoint indices based on the desired processing
% % length and degree of overlap
% for i = 1:length(break_idx)-1
%     
%     % Determine lat/lon and distance along the ith section radar data
%     lat_i = lat(break_idx(i)+1:break_idx(i+1));
%     lon_i = lon(break_idx(i)+1:break_idx(i+1));
%     dist_i = pathdist(lat_i, lon_i);
%     
%     % Initialize while loop for current radar segment
%     search = true;
%     j_start = break_idx(i)+1;
%     
%     % If the entire radar segment is smaller than the desired minimum
%     % length, add breakpoints to process a radargram over the entire
%     % segment
%     if dist_i(end) <= length_min
%         search = false;
%         j_end = length(dist_i);
%         breaks_new{j} = [j_start j_end+break_idx(i)];
%         j = j+1;
%     end
%     
%     % Iteratively find breakpoint indices with desired radargram length and
%     % overlap within ith continuous data segment
%     while search == true
%         
%         % Calculate distances relative to current starting breakpoint
%         dist_j = dist_i - dist_i(j_start-break_idx(i));
%         
%         % Find ending breakpoint index (relative to distances along current
%         % data segment)
%         j_end = find(dist_j >= length_min, 1);
%         
%         % Check if only minor additional data (less than the overlap
%         % distance) remains in current data segment
%         if dist_j(end)-dist_j(j_end) <= overlap
%             
%             % If only minor data remains, add remainder to current
%             % radargram breakpoints and end search for additional
%             % breakpoints within ith data segment
%             j_end = length(dist_j);
%             search = false;
%         end
%         
%         % Check if remaining data is shorter than minimum desired length
%         if isempty(j_end)
%             
%             % If remainder is too short, set starting breakpoint such that
%             % desired length is achieved
%             j_start = find(dist_j>=dist_j(end)-length_min, 1) + break_idx(i);
%             
%             % Set ending breakpoint to end of current data segment, and end
%             % search for addtional breakpoints
%             j_end = length(dist_i);
%             search = false;
%         end
%         
%         % Export breakpoints to preallocated array (use absolute indices,
%         % not relative to current data segment)
%         breaks_new{j} = [j_start j_end+break_idx(i)];
%         
%         % Set starting breakpoint position for next while iteration based
%         % on the ending breakpoint and overlap distance
%         j_start = find(dist_j-dist_j(j_end)+overlap>=0, 1) + break_idx(i);
%         j = j+1;
%     end
% end
% 
% %%
% 
% % Preallocate arrays for file indices and trace position indices to include
% % in each continuous radargram
% files_i = zeros(length(breaks_new), 2);
% position_i = zeros(length(breaks_new), 2);
% 
% % Loop through each set of breakpoints to determine which files and traces
% % to include in each continuous radargram
% for i = 1:length(breaks_new)
%     
%     % Find the file index at start of ith continuous radargram using the
%     % previously calculated breakpoint indices
%     files_i(i,1) = find(file_end_idx>=breaks_new{i}(1), 1);
%     
%     % Find the file index at the end of ith continuous radargram using the
%     % previously calculated breakpoint indices
%     files_i(i,2) = find(file_end_idx>=breaks_new{i}(2), 1);
%     
%     % Find the trace indices within the first and last file corresponding
%     % to the intra-file breakpoints of the ith continuous radargram
%     position_i(i,1) = breaks_new{i}(1) - ...
%         (file_end_idx(files_i(i,1)) - file_end_idx(1));
%     position_i(i,2) = breaks_new{i}(2) - (file_end_idx(files_i(i,2)-1));
%     
% %     check1_idx = file_end_idx(files_i(i,1)) - file_end_idx(1) + position_i(i,1)
% %     d(check1_idx)
% %     check2_idx = file_end_idx(files_i(i,2)-1) + position_i(i,2)
% %     d(check2_idx)
% end
% 
% %%
% 
% for i = 1:size(files_i,1)
%     
%     switch source
%         case 'SEAT'
%             % Initialize for loop with first file in directory
%             data_struct = import_radar(fullfile(files(files_i(i,1)).folder, ...
%                 files(files_i(i,1)).name), position_i(i,1));
%         case 'OIB'
%             % Initialize for loop with first file in directory
%             OIB_data = OIB_convert(fullfile(files(files_i(i,1)).folder, ...
%                 files(files_i(i,1)).name));
%             data_struct = import_radar(OIB_data, position_i(i,1));
%     end
%     
%     % Convert first imported file to cell array
%     data2cells = struct2cell(data_struct);
%     
%     % Get field names in data structure
%     fld_names = fieldnames(data_struct);
%     
%     % List of desired field names to search for and include in radar
%     % processing
%     fld_wanted = {'collect_date', 'lat', 'lon', 'elev', 'Easting', 'Northing',...
%         'dist', 'data_out', 'arr_layers', 'time_trace'};
%     
%     % Preallocate logical array for field names to include
%     fld_include = logical(zeros(length(fld_names),1));
%     
%     % Find the indices of field names in the data structure which match the
%     % desired field names to include in processing
%     for j = 1:length(fld_include)
%         fld_include(j) = any(cellfun(@(x) strcmpi(fld_names{j}, x), fld_wanted));
%     end
%     data_cells = data2cells(fld_include);
%     
%     % Concatenate adjacent radar files within the current radargram
%     for j = files_i(i,1)+1:files_i(i,2)-1
%         switch source
%             case 'SEAT'
%                 % Import and clean each radar file
%                 data_full = struct2cell(import_radar(fullfile(files(j).folder, ...
%                     files(j).name)));
%             case 'OIB'
%                 % Import and clean each radar file
%                 OIB_data = OIB_convert(fullfile(files(j).folder, ...
%                     files(j).name));
%                 data_full = struct2cell(import_radar(OIB_data));
%         end
%         
%         % Select only the desired fields within current cells
%         data_j = data_full(fld_include);
%         
%         % Add data from cells to the overall radargram data cells
%         data_cells = cellfun(@horzcat, data_cells, data_j, 'UniformOutput', 0);
%     end
%     
%     
%     switch source
%         case 'SEAT'
%             % Import and clean each radar file
%             data_full = struct2cell(import_radar(fullfile(...
%                 files(files_i(i,2)).folder, files(files_i(i,2)).name), ...
%                 [], position_i(i,2)));
%         case 'OIB'
%             % Import and clean each radar file
%             OIB_data = OIB_convert(fullfile(files(files_i(i,2)).folder, ...
%                 files(files_i(i,2)).name));
%             data_full = struct2cell(import_radar(OIB_data, ...
%                 [], position_i(i,2)));
%     end
%     
%     % Select only the desired fields within the final cell array in the ith
%     % radargram
%     data_j = data_full(fld_include);
%     
%     % Add data from the final cell array to the ith radargram data cells
%     data_cells = cellfun(@horzcat, data_cells, data_j, 'UniformOutput', 0);
%     
%     % Average values for collect_date for ith continuous radargram
%     date_idx = cellfun(@(x) strcmpi('collect_date', x), fld_names(fld_include));
%     data_cells{date_idx} = median(data_cells{date_idx}, 2);
%     
%     % Average values for time_trace for ith continuous radargram
%     time_idx = cellfun(@(x) strcmpi('time_trace', x), fld_names(fld_include));
%     data_cells{time_idx} = median(data_cells{time_idx}, 2);
%     
% %     % Average values for TWTT for ith continuous radargram (should only
% %     % exist at this point for OIB data)
% %     TWTT_idx = cellfun(@(x) strcmpi('TWTT', x), fld_names(fld_include));
% %     try
% %         data_cells{TWTT_idx} = median(data_cells{TWTT_idx}, 2);
% %     catch
% %         disp('Warning: missing TWTT data; data should be SEAT radar')
% %     end
%     
%     % Convert concatenated cell arrays back to structure with desired
%     % fieldnames
%     radar_tmp = cell2struct(data_cells, fld_names(fld_include), 1);
%     
%     % Recalculate continuous path distance along the concatenated ith
%     % radargram
%     radar_tmp.dist = pathdist(radar_tmp.lat, radar_tmp.lon);
%     
%     tic
%     [radar_tmp] = radar_RT(radar_tmp, cores, Ndraw);
%     [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
%     toc
%     
%     
% end
