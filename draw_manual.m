% Script to load selected radargram and manually draw annual layers, and
% then save resulting layer positions for later use

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
        computer = 'laptop';
        %         computer = input('Current PC: ');
        switch computer
            case 'work'
                data_path = 'E:/Research/Antarctica/Data/';
                addon_path = 'C:/Users/u1046484/Documents/MATLAB/Addons/';
                
            case 'laptop'
                %                 data_path = 'C:/Users/durba/Documents/Research/Antarctica/Data/';
                data_path = 'F:/Research/Antarctica/Data/';
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

%% 

% Name of SEAT core site to generate training data/perform regression
name = 'SEAT10_4';

% Load relevant radar data (previously generated using the above section)
radar = load(fullfile(data_path, 'IceBridge/manual_layers', name, ...
    strcat('layers_', name, '.mat')));

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
    hi = drawpolyline();
    
    if isvalid(hi)
        
        % Find the range of the manually picked layer
        col = (1:length(radar.Easting))';
%         col = (max([1 round(min(hi.Position(:,1)))]):...
%             min([round(max(hi.Position(:,1))) length(radar.Easting)]))';
        
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

man_layers = man_layers(~cellfun(@isempty,man_layers));