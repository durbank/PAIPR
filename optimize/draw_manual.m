% Script to load selected radargram and manually draw annual layers, and
% then save resulting layer positions for later use

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:/Research/Antarctica/Data/';
        addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
    case false
        data_path = '/media/durbank/WARP/Research/Antarctica/Data/';
        addon_path = '/home/durbank/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_struct = dir(fullfile(addon_path, 'AntarcticMappingTools_*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add export_fig to path
addon_struct = dir(fullfile(addon_path, 'altmany-export_fig*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))
% Add CReSIS OIB MATLAB reader functions to path
addon_struct = dir(fullfile(addon_path, 'cresis-L1B-matlab-readers*'));
addpath(genpath(fullfile(addon_struct.folder, addon_struct.name)))

% Add PAIPR-core functions to path
% parent = cd;
% parent = fullfile(parent,'..', 'PAIPR-core');
PAIPR_path = fullfile(cd,'..', 'PAIPR-core');
addpath(genpath(PAIPR_path))

%% Process raw OIB echograms for PAIPR results and save output

% % Define number of Monte Carlo simulations to perform
% Ndraw = 100;
% 
% % Load core data from file (data used was previously generated using
% % import_cores.m)
% core_file = fullfile(data_path, 'Ice-cores/SEAT_cores/SEAT_cores.mat');
% cores = load(core_file);
% 
% % Select radar directory for processing
% [input_dir] = uigetdir(data_path,...
%     "Select directory containing raw echograms to process");
% 
% % Select output directory in which to save processed echogram and manual
% % layers
% [output_dir] = uigetdir(input_dir, ...
%     "Select directory to output processed echogram");
% 
% % Process OIB echogram with PAIPR
% [radar] = PAIPR_draw(input_dir, cores, Ndraw);
% 
% % Save processed radar structure for future use
% output_path = fullfile(output_dir, "PAIPR_out.mat");
% save(output_path, '-struct', 'radar', '-v7.3')

%% Load previously processed PAIPR echogram to manually trace layers

% Load relevant radar data (previously generated using the above section)
[r_file, r_path] = uigetfile(data_path,...
    "Select radar file to use in manual picking");
radar_file = fullfile(r_path, r_file);
radar = load(radar_file);

% guides = load(fullfile(data_path, 'IceBridge/manual_layers', name, ...
%     strcat('manual_', name, '.mat')));
% guides = guides.man_all;

%%
% Figure for manually tracing visible annual layers

% Preallocate cell array for position subscripts of manual layers
man_layers = cell(1,250);
i = 1;
draw = true;
while draw==true
    
    % Draw position of manual layers in radargram (layers should be continuous
    % across the entire radargram)
    f_draw = figure;
    imagesc(radar.data_smooth, [-2 2])
    hold on
%     for j = 1:length(guides)
%         plot(guides{j}(:,1), guides{j}(:,2), 'm')
%     end
    if i >= 2
        for k = 1:i-1
            plot(man_layers{k}(:,1), man_layers{k}(:,2), 'r')
        end
    end
    hi = drawpolyline();
    
    if isvalid(hi)
        
        % Find the range of the manually picked layer
%         col = (1:length(radar.Easting))';
        col = (max([1 round(min(hi.Position(:,1)))]):...
            min([round(max(hi.Position(:,1))) length(radar.Easting)]))';
        
        % Linearly interpolate layer row positions to the full range of the
        % manually picked layer
        row = interp1(hi.Position(:,1), hi.Position(:,2), col, ...
            'linear', 'extrap');
        
        % Export layer position subscripts to preallocated cell array
        man_layers{i} = [col row];
        i = i+1;
        
        close(f_draw)
        clear hi
    else
        draw = false;
    end
    
end

% Remove empty cells
man_layers = man_layers(~cellfun(@isempty,man_layers));

% Keep only traced layer data inside the boundaries of the echogram
keep_idx = cellfun(@(x) round(x(:,2))<=size(radar.data_smooth,1), ...
    man_layers, 'UniformOutput', false);
man_layers = cellfun(@(x,y) x(y,:), man_layers, keep_idx, 'UniformOutput', false);


%% Save manual layer output to disk for later use

% Save processed radar structure for future use
output_path = fullfile(output_dir, "manual_layers.mat");
save(output_path, '-struct', 'man_layers.man_all', '-v7.3')
