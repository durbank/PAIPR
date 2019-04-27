% Function to format and prepare for export PAIPR-processed data to various
% different forms



% GUI to specify input data directory
input_dir = uigetdir('../', ...
    'Select directory from which to input PAIPR-processed data');

% Find file names for previously processed PAIPR data
wild = '*.mat';
files = dir(fullfile(input_dir, wild));

% Preallocate arrays of sufficient size for data
Easting = zeros(1, length(files)*2*(50*1000/25));
Northing = Easting;
SMB = cell(1, length(files)*2*(50*1000/25));
year = SMB;

for i = 1:length(files)
    
    % Load relevent data from current data file
    load(fullfile(files(i).folder, files(i).name), 'Easting');
    load(fullfile(files(i).folder, files(i).name), 'Northing');
    load(fullfile(files(i).folder, files(i).name), 'SMB');
    load(fullfile(files(i).folder, files(i).name), 'SMB_yr');
    
    % Find position of last data entered into preallocated arrays
    next_idx = sum(~cellfun(@isempty, SMB)) + 1;
    
    % Fill current iteration data into preallocated arrays
    Easting(next_idx:next_idx+length(Easting)-1) = Easting;
    Northing(next_idx:next_idx+length(Northing)-1) = Northing;
    SMB(next_idx:next_idx+length(SMB)-1) = SMB;
    year(next_idx:next_idx+length(SMB_yr)-1) = SMB_yr;
    
end

% Find and remove empty SEAT indices (usually from excess preallocation)
keep_idx = find(~cellfun(@isempty, SMB));
Easting = Easting(keep_idx);
Northing = Northing(keep_idx);
SMB = SMB(keep_idx);
year = year(keep_idx);

% Attempt to additionally load elevation data, if available
try
    elev = false(1, length(files)*2*(50*1000/25));
    for i=1:length(files)
        load(fullfile(files(i).folder, files(i).name), 'elev');
        next_idx = sum(logical(elev)) + 1;
        elev(next_idx:next_idx+length(elev)-1) = elev;
    end
    keep_idx = find(logical(elev));
    elev = elev(keep_idx);
catch
    disp('Flag: Missing elevation data')
end