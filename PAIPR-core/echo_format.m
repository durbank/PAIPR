% Function to import and process an arbitrary number of radar echogram
% files within a single directory
% File formats may either be in NASA Operation IceBridge Snow radar
% standard format (.nc) or in SEAT Ku-band radar echogram standard format
% (.mat)

function [radar_ALL] = echo_format(radar_dir)

% List all files matching 'wild' within radar directory
wild = '*.mat';
files = dir(fullfile(radar_dir, wild));

format = 'MATLAB';
source = 'SEAT';

% If directory has no '.mat' files, additionally check for '.nc' files
if isempty(files)
    wild = '*.nc';
    files = dir(fullfile(radar_dir, wild));
    format = 'NetCDF';
    source = 'OIB';
end

%%

% Preallocate arrays for continuous lat/lon/time positions for all data in
% directory
lat = [];
lon = [];
time = [];

% Preallocate array for the radargram length (in data bins) for each
% data file in directory
file_length = zeros(1, length(files));

switch format
    case 'MATLAB'
        
        % Iteratively load lat/lon positions for each data file, and add
        % corresponding values to arrays
        for i = 1:length(files)
            data_i = load(fullfile(files(i).folder, files(i).name), ...
                'lat', 'lon');
            file_length(i) = length(data_i.lat);
            lat = [lat data_i.lat];
            lon = [lon data_i.lon];
%             time = [time data_i.time_gps];
        end
        
    case 'NetCDF'
        
        % Iteratively load lat/lon positions for each data file, and add
        % corresponding values to arrays
        for i = 1:length(files)
            lat_i = ncread(fullfile(files(i).folder, files(i).name), 'lat');
            lon_i = ncread(fullfile(files(i).folder, files(i).name), 'lon');
            time_i = ncread(fullfile(files(i).folder, files(i).name), 'time');
            file_length(i) = length(lat_i);
            lat = [lat lat_i'];
            lon = [lon lon_i'];
            time = [time time_i'];
        end
        
end

% Arrange lat/lon data in ascending order of collection time
if ~isempty(time)
    [~,sort_idx] = sort(time, 'ascend');
    lat = lat(sort_idx);
    lon = lon(sort_idx);
end

% Calculate the cummulative radargram length (in data bins) for the data
% files in directory
file_end_idx = [0 cumsum(file_length)];
% file_end_idx = cumsum(file_length);

% Replace data without valid location values (outside Antarctic Circle)
% with nearest preceding valid location (these missing data will later be
% removed in processing)
invld_idx = lat>=-65 | lat<-90;
lat(invld_idx) = NaN;
lon(invld_idx) = NaN;
lat = fillmissing(lat, 'previous');
lon = fillmissing(lon, 'previous');

% Calucate distance along traverse (in meters)
d = pathdist(lat, lon);

% Find indices to break radargrams based on absence of data across an
% extended distance (greater than 500 m)
break_idx = [0 find(diff(d)>500) length(lat)];

% Set the minimum length needed for radargram processing and radargram
% overlap interval (in meters)
length_min = 30000;
overlap = 10000;

% Preallocate and initialize variables for breakpoint calculation
breaks_new = cell(1, length(break_idx)-1);
j = 1;

% For each continuous section of radargram (no significant breaks),
% determine additional breakpoint indices based on the desired processing
% length and degree of overlap
for i = 1:length(break_idx)-1
    
    % Determine lat/lon and distance along the ith section radar data
    lat_i = lat(break_idx(i)+1:break_idx(i+1));
    lon_i = lon(break_idx(i)+1:break_idx(i+1));
    dist_i = pathdist(lat_i, lon_i);
    
    % Initialize while loop for current radar segment
    search = true;
    j_start = break_idx(i)+1;
    
    % If the entire radar segment is smaller than the desired minimum
    % length, add breakpoints to process a radargram over the entire
    % segment
    if dist_i(end) <= length_min
        search = false;
        j_end = length(dist_i);
        breaks_new{j} = [j_start j_end+break_idx(i)];
        j = j+1;
    end
    
    % Iteratively find breakpoint indices with desired radargram length and
    % overlap within ith continuous data segment
    while search == true
        
        % Calculate distances relative to current starting breakpoint
        dist_j = dist_i - dist_i(j_start-break_idx(i));
        
        % Find ending breakpoint index (relative to distances along current
        % data segment)
        j_end = find(dist_j >= length_min, 1);
        
        % Check if only minor additional data (less than the overlap
        % distance) remains in current data segment
        if dist_j(end)-dist_j(j_end) <= overlap
            
            % If only minor data remains, add remainder to current
            % radargram breakpoints and end search for additional
            % breakpoints within ith data segment
            j_end = length(dist_j);
            search = false;
        end
        
        % Check if remaining data is shorter than minimum desired length
        if isempty(j_end)
            
            % If remainder is too short, set starting breakpoint such that
            % desired length is achieved
            j_start = find(dist_j>=dist_j(end)-length_min, 1) + break_idx(i);
            
            % Set ending breakpoint to end of current data segment, and end
            % search for addtional breakpoints
            j_end = length(dist_i);
            search = false;
        end
        
        % Export breakpoints to preallocated array (use absolute indices,
        % not relative to current data segment)
        breaks_new{j} = [j_start j_end+break_idx(i)];
        
        % Set starting breakpoint position for next while iteration based
        % on the ending breakpoint and overlap distance
        j_start = find(dist_j-dist_j(j_end)+overlap>=0, 1) + break_idx(i);
        j = j+1;
    end
end

%%

% Preallocate arrays for file indices and trace position indices to include
% in each continuous radargram
files_i = zeros(length(breaks_new), 2);
position_i = zeros(length(breaks_new), 2);

% Loop through each set of breakpoints to determine which files and traces
% to include in each continuous radargram
for i = 1:length(breaks_new)
    
    % Find the file index at start of ith continuous radargram using the
    % previously calculated breakpoint indices
    files_i(i,1) = find(file_end_idx>=breaks_new{i}(1), 1) - 1;
    
    % Find the file index at the end of ith continuous radargram using the
    % previously calculated breakpoint indices
    files_i(i,2) = find(file_end_idx>=breaks_new{i}(2), 1) - 1;
    
    % Find the trace indices within the first and last file corresponding
    % to the intra-file breakpoints of the ith continuous radargram
    position_i(i,1) = breaks_new{i}(1) - file_end_idx(files_i(i,1));
    position_i(i,2) = breaks_new{i}(2) - file_end_idx(files_i(i,2));
end

%%

radar_ALL = struct;

for i = 1:size(files_i,1)
    
    switch source
        case 'SEAT'
            % Initialize for loop with first file in directory
            data_struct = orderfields(import_radar(fullfile(...
                files(files_i(i,1)).folder, files(files_i(i,1)).name), ...
                position_i(i,1)));
        case 'OIB'
            % Initialize for loop with first file in directory
            OIB_data = OIB_convert(fullfile(files(files_i(i,1)).folder, ...
                files(files_i(i,1)).name));
            data_struct = orderfields(import_radar(OIB_data, position_i(i,1)));
    end
    
    % Convert first imported file to cell array
    data2cells = struct2cell(data_struct);
    
    % Get field names in data structure
    fld_names = fieldnames(data_struct);
    
    % List of desired field names to search for and include in radar
    % processing
    fld_wanted = {'collect_date', 'lat', 'lon', 'elev', 'Easting', ...
        'Northing', 'time_gps', 'dist', 'data_out', 'time_trace'};
%     fld_wanted = {'collect_date', 'lat', 'lon', 'elev', 'Easting', 'Northing',...
%         'dist', 'data_out', 'arr_layers', 'time_trace'};

    % Preallocate logical array for field names to include
    fld_include = zeros(length(fld_names),1);
    
    % Find the indices of field names in the data structure which match the
    % desired field names to include in processing
    for j = 1:length(fld_include)
        fld_include(j) = any(cellfun(@(x) strcmpi(fld_names{j}, x), fld_wanted));
    end
    fld_include = logical(fld_include);
    data_cells = data2cells(fld_include);
    
    % Concatenate adjacent radar files within the current radargram
    for j = files_i(i,1)+1:files_i(i,2)-1
        switch source
            case 'SEAT'
                % Import and clean each radar file
                data_full = struct2cell(orderfields(import_radar(...
                    fullfile(files(j).folder, files(j).name))));
            case 'OIB'
                % Import and clean each radar file
                OIB_data = OIB_convert(fullfile(files(j).folder, ...
                    files(j).name));
                data_full = struct2cell(orderfields(import_radar(OIB_data)));
        end
        
        % Select only the desired fields within current cells
        data_j = data_full(fld_include);
        
        % Add data from cells to the overall radargram data cells
        data_cells = cellfun(@horzcat, data_cells, data_j, 'UniformOutput', 0);
    end
    
    
    switch source
        case 'SEAT'
            % Import and clean each radar file
            data_full = struct2cell(orderfields(import_radar(fullfile(...
                files(files_i(i,2)).folder, files(files_i(i,2)).name), ...
                [], position_i(i,2))));
        case 'OIB'
            % Import and clean each radar file
            OIB_data = OIB_convert(fullfile(files(files_i(i,2)).folder, ...
                files(files_i(i,2)).name));
            data_full = struct2cell(orderfields(import_radar(OIB_data, ...
                [], position_i(i,2))));
    end
    
    % Select only the desired fields within the final cell array in the ith
    % radargram
    data_j = data_full(fld_include);
    
    % Add data from the final cell array to the ith radargram data cells
    data_cells = cellfun(@horzcat, data_cells, data_j, 'UniformOutput', 0);
    
    % Average values for collect_date for ith continuous radargram
    date_idx = cellfun(@(x) strcmpi('collect_date', x), fld_names(fld_include));
    data_cells{date_idx} = median(data_cells{date_idx}, 2);
    
    % Average values for time_trace for ith continuous radargram
    time_idx = cellfun(@(x) strcmpi('time_trace', x), fld_names(fld_include));
    data_cells{time_idx} = median(data_cells{time_idx}, 2);
    
%     % Average values for TWTT for ith continuous radargram (should only
%     % exist at this point for OIB data)
%     TWTT_idx = cellfun(@(x) strcmpi('TWTT', x), fld_names(fld_include));
%     try
%         data_cells{TWTT_idx} = median(data_cells{TWTT_idx}, 2);
%     catch
%         disp('Warning: missing TWTT data; data should be SEAT radar')
%     end
    
    % Convert concatenated cell arrays back to structure with desired
    % fieldnames
    radar_tmp = cell2struct(data_cells, fld_names(fld_include), 1);
    
    % Recalculate continuous path distance along the concatenated ith
    % radargram
    radar_tmp.dist = pathdist(radar_tmp.lat, radar_tmp.lon);
    
    radar_ALL(i).segment = radar_tmp;
end

end