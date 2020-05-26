% Function to process OIB Snow echograms with PAIPR on CHPC using SLURM

% DEPENDENCIES: Requires 'PAIPR-core' functions, 'cresis-L1B-matlab-reader'
%               functions, and 'AntarcticMappingTools' functions to be
%               added to the path previously

function [success_codes] = process_PAIPR(...
    DATADIR, RHO_PATH, OUTDIR, Ndraw, NameValueArgs)

arguments
    DATADIR %string
    RHO_PATH %string
    OUTDIR %string
    Ndraw %int16
    NameValueArgs.VerboseOutput = false
end

% % Define number of Monte Carlo simulations to perform
% Ndraw = 100;

% Directory containing raw OIB echograms
radar_dir = DATADIR;

% % Select modeled density subset .mat file to use, and load to workspace
% rho_file = load(RHO_PATH);
% rho_subset = rho_file.rho_subset;

echo_out = fullfile(OUTDIR, 'echo');
if NameValueArgs.VerboseOutput == true
    % Check for existence of echogram output dir and create as necessary
    if 7~=exist(echo_out, 'dir')
        mkdir(echo_out)
    end
end

% Check for existence of .csv output directory and create as necessary
gamma_out = fullfile(OUTDIR, 'gamma');
if 7~=exist(gamma_out, 'dir')
    mkdir(gamma_out)
end

%%

% Import density profiles and format as nested table
rho_table = import_rho(RHO_PATH);

% Get grouping indices of echograms
[files, start_idx, end_idx] = OIB_chunk(radar_dir);


% Define overlap distance and the final horizontal resolution of the output
% data
horz_res = 25;
overlap_switch = false;
overlap = 10000;


echo_0 = orderfields(import_radar(...
    OIB_convert(fullfile(files(1).folder, files(1).name))));
fields = fieldnames(echo_0);
for i = 1:length(fields)
    echo_0.(fields{i}) = [];
end

% Preallocate array for success codes
success_codes = zeros(1, length(end_idx));

parfor i=1:length(end_idx)
    
    try
        echo_i = echo_0;
        for j = start_idx(i):end_idx(i)
            
            OIB_j = OIB_convert(fullfile(files(j).folder, files(j).name));
            struct_j = orderfields(import_radar(OIB_j));
            
            echo_i = cell2struct(cellfun(@horzcat, struct2cell(echo_i), ...
                struct2cell(struct_j), 'uni', 0), fieldnames(echo_i), 1);
            
            
        end
        echo_i.dist = pathdistps(echo_i.lat, echo_i.lon);
        
        
        
        [radar_tmp] = radar_stack(echo_i, horz_res);
        
        
        % Load modeled depth-density data from stats model output at 
        % specified Easting/Northing coordinates
        [rho_data] = load_rho(...
            rho_table, radar_tmp.Easting, radar_tmp.Northing);
        
        % Convert to depth
        [radar_tmp] = radar_depth(radar_tmp, rho_data);
        
        % Calculate radar age-depth profile distributions (includes 
        % processing signal-noise, radon transforms, layer tracing, 
        % likelihood assignments, and age calculations)
        [radar_tmp] = calc_layers(radar_tmp, 'stream');
        
        % Calculate age-depth profiles
        r = -6.723;
        k = 3.0;
        [radar_tmp] = radar_age(radar_tmp, r, k, Ndraw);
        
        % Calculate radar annual SMB
        [radar_tmp] = calc_SWE(radar_tmp, rho_data, Ndraw);
        
        % Perform QC check on echogram image
        [QC_med, QC_val, QC_flag, depth_idx, yr_cutoff] = QC_check(...
            radar_tmp, 0.50, 0.10);
        
        switch overlap_switch
            case true
                % Clip radar data structure variables based on the desired 
                % radargram %overlap, and combine desired clipped variables 
                % into new data structure
                clip = round(0.5*overlap/horz_res);
                radar = struct(...
                    'collect_time',radar_tmp.collect_time(clip:end-clip),...
                    'Easting', radar_tmp.Easting(clip:end-clip),...
                    'Northing', radar_tmp.Northing(clip:end-clip), ...
                    'dist', radar_tmp.dist(clip:end-clip), ...
                    'depth', radar_tmp.depth, ...
                    'data_smooth', radar_tmp.data_smooth(:,clip:end-clip),...
                    'IM_grad', radar_tmp.IM_grad(:,clip:end-clip), ...
                    'IM_QC', radar_tmp.IM_QC(:,clip:end-clip), ...
                    'QC_med', QC_med, 'QC_val', QC_val, ...
                    'QC_flag', QC_flag', 'QC_depth_idx', depth_idx,...
                    'QC_yr', yr_cutoff, ...
                    'DB', radar_tmp.DB(:,clip:end-clip), ...
                    'likelihood', radar_tmp.likelihood(:,clip:end-clip), ...
                    'ages', radar_tmp.ages(:,clip:end-clip,:));
                
                
                if isfield(radar_tmp, 'groups')
                    radar.groups = radar_tmp.groups(:,clip:end-clip);
                    
                    % Find layer member indices based on new clipped record
                    layers = cell(1, max(radar.groups(:)));
                    for j = 1:length(layers)
                        layers{j} = find(radar.groups == j);
                    end
                    radar.layers = layers(~cellfun(@isempty, layers));
                    
                end
                
                % Check for existing elevation data, and add to output 
                % struct if available
                if isfield(radar_tmp, 'elev')
                    radar.elev = radar_tmp.elev(clip:end-clip);
                else
                    radar.elev = nan(1, length(radar.Easting));
                end
                
                % Redefine radargram distances based on clipped data
                radar.dist = radar.dist - radar.dist(1);
                
                % Clip SMB cell array data and add to output struct
                radar.SMB_yr =  radar_tmp.SMB_yr(clip:end-clip);
                radar.SMB = radar_tmp.SMB(clip:end-clip);
                
            case false
                
                radar = radar_tmp;
                radar.QC_med = QC_med;
                radar.QC_val = QC_val;
                radar.QC_flag = QC_flag;
                radar.QC_depth_idx = depth_idx;
                radar.QC_yr = yr_cutoff;
        end
        
        % Generate file names and paths under which to save data
        filename = sprintf('%s%d','radar_out',i);
        csv_output = fullfile(gamma_out, strcat(filename, '.csv'));
        
        switch NameValueArgs.VerboseOutput
            case true
                mat_output = fullfile(echo_out, strcat(filename, '.mat'));
                % Save output structures to disk
                [success_codes(i)] = parsave(...
                    radar, csv_output, mat_output);
            case false
                % Save output structures to disk
                [success_codes(i)] = parsave(radar, csv_output);
        end
        
    catch ME
        fprintf(1, ['An error occurred while processing data in the '...
        'following directory: \n %s \n'], files(i).folder);
        fprintf(1, "Affected echograms span %s to %s \n", ...
            files(start_idx(i)).name, files(end_idx(i)).name)
        fn_missing = sprintf('%s%d', 'radar_out', i);
        fprintf(1, "Resultant missing files should be named %s \n", ...
            fn_missing);
        fprintf(1,"The error occurred in the following location: %s \n",...
            ME.stack)
        fprintf(1, "The error thrown reads:\n%s \n", ME.message);
    end
    
    
end

end