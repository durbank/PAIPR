

% Start parellel pool
poolobj=parpool('local', 16);

%%
% Add PAIPR-core functions to path
addon_struct = dir(fullfile('src/', 'PAIPR-core*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add cresis-matlab-reader functions to path
addon_struct = dir(fullfile('src/Dependencies/', ...
    'cresis-L1B-matlab-reader*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile('src/Dependencies/', ...
    'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

%     % Add PAIPR scripts to path
%     addpath('src/scripts/')

%%

% Define input/output file locations
%rho_file = 'rho_data.csv';
%     rho_files = dir(fullfile('rho_data', '*.csv'));
rho_dir = 'rho_data';
out_dir = 'Outputs/';

% Define directories for raw echogram inputs
data_dirs = dir('Data');
data_dirs = data_dirs([data_dirs(:).isdir]);
datadir_list = data_dirs(~ismember({data_dirs(:).name},{'.','..'}));


for i=1:length(datadir_list)
    
    try
        
        % Get raw echogram directory for current iteration
        echo_dir = datadir_list(i);
        
        % Create output directory for current iteration
        out_i = fullfile(out_dir, echo_dir.name);
        mkdir(out_i);
        
        % Get rho file for current iteration
        rho_file = fullfile(rho_dir, ...
            strcat('rho_', echo_dir.name, '.csv'));
        
        if ~isfile(rho_file)
            rho_file = fullfile('src/data-defaults', 'density-default.csv');
        end
        
        % Run PAIPR functions
        NSIM = 250;
        [success_codes] = process_PAIPR(...
            fullfile(echo_dir.folder, echo_dir.name), ...
            rho_file, out_i, NSIM, ...
            'LayerDepths', true);
        %         [success_codes] = process_PAIPR(...
        %             fullfile(echo_dir.folder, echo_dir.name), ...
        %             rho_file, out_i, NSIM, 'VerboseOutput', true);
        
        % Save sucess codes as .csv table
        T = table(1:length(success_codes), success_codes, ...
            'VariableNames', {'Iteration', 'Exit_code'});
        writetable(T, fullfile(out_i, 'exit_codes.csv'))
        
    catch ME
        disp("An error occurred in processing");
        disp("The error occurred in the following locations: ")
        for m=1:length(ME.stack)
            fprintf(1, '>> Function = %s: Line %u\n', ...
                ME.stack(m).name, ME.stack(m).line)
        end
        fprintf(1, "The error thrown reads:\n%s\n", ME.message);
        disp('--------------------------------------------------')    
    end
    
end

%%

% If required, close parellel pool
delete(poolobj)

% Exit MATLAB server
exit
