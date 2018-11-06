% Script compare results from auto-IceBridge to cores and manual IceBridge

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

% name = 'SEAT10_5';
% input_dir = fullfile(data_path, 'IceBridge/manual_layers', name);
% radar_ALL = radar_format(fullfile(input_dir, 'raw_data/'));
% radar = radar_ALL(1).segment;
% overlap = 10000;
% horz_res = 25;
% 
% [radar_tmp] = radar_RT(radar, cores, Ndraw);
% [radar_tmp] = calc_SWE(radar_tmp, Ndraw);
% % 
% clip = round(0.5*overlap/horz_res);
% 
% fld_nm = fieldnames(radar_tmp);
% fld_want = {'collect_date', 'Easting', 'Northing', 'dist', 'depth', ...
%     'rho_coeff', 'rho_var', 'data_smooth', 'peaks', 'groups', 'ages', ...
%     'SMB_yr', 'SMB'};
% 
% radar = struct('collect_date', radar_tmp.collect_date, ...
%     'Easting', radar_tmp.Easting(clip:end-clip),...
%     'Northing', radar_tmp.Northing(clip:end-clip), ...
%     'dist', radar_tmp.dist(clip:end-clip), 'depth', radar_tmp.depth, ...
%     'rho_coeff', radar_tmp.rho_coeff(:,clip:end-clip), ...
%     'rho_var', radar_tmp.rho_var(:,clip:end-clip), ...
%     'data_smooth', radar_tmp.data_smooth(:,clip:end-clip),...
%     'peaks', radar_tmp.peaks(:,clip:end-clip), ...
%     'groups', radar_tmp.groups(:,clip:end-clip),...
%     'likelihood', radar_tmp.likelihood(:,clip:end-clip), ...
%     'ages', radar_tmp.ages(:,clip:end-clip,:));
% radar.dist = radar.dist - radar.dist(1);
% if isfield(radar_tmp, 'elev')
%     radar.elev = radar_tmp.elev(clip:end-clip);
% end
% radar.SMB_yr =  radar_tmp.SMB_yr(clip:end-clip);
% radar.SMB = radar_tmp.SMB(clip:end-clip);
% % 
% fn = strcat('layers_', name, '_opt.mat');
% output_path = fullfile(input_dir, fn);
% save(output_path, '-struct', 'radar', '-v7.3')


%%

% Name of SEAT core site to generate training data/perform regression
name = 'SEAT10_6';

% Load relevant radar data (previously generated using the above section)
radar = load(fullfile(data_path, 'IceBridge/manual_layers', name, ...
    strcat('layers_', name, '_opt.mat')));

% Load manually traced layers for current SEAT core site (generated using
% 'draw_manual.m')
tmp = load(fullfile(data_path, 'IceBridge/manual_layers', name, ...
    strcat('manual_', name, '.mat')));
man_layers = tmp.man_all;

keep_idx = cellfun(@(x) round(x(:,2))<=size(radar.data_smooth,1), ...
    man_layers, 'UniformOutput', false);
man_layers = cellfun(@(x,y) x(y,:), man_layers, keep_idx, 'UniformOutput', false);

%%

man_mat = zeros(size(radar.data_smooth));

for i = 1:length(man_layers)
    
    layer_idx = sub2ind(size(man_mat), round(man_layers{i}(:,2)), man_layers{i}(:,1));
    man_mat(layer_idx) = i;
end




man_log = logical(man_mat);
ages_man = zeros(size(radar.data_smooth));
age_top = radar.collect_date;
yr_pick1 = ceil(radar.collect_date - 1);

for i = 1:size(ages_man,2)
    depths_i = [0; radar.depth(man_log(:,i))];
    yrs_i = ([age_top yr_pick1:-1:yr_pick1-length(depths_i)+2])';
    ages_man(:,i) = interp1(depths_i, yrs_i, radar.depth, 'linear', 'extrap');
end


%%


figure
hold on
plot(radar.depth, mean(ages_man,2), 'k', 'LineWidth', 2)
plot(radar.depth, mean(ages_man,2) + std(ages_man, [], 2), 'k--')
plot(radar.depth, mean(ages_man,2) - std(ages_man, [], 2), 'k--')
plot(radar.depth, mean(mean(radar.ages, 3), 2), 'm', 'LineWidth', 2)
plot(radar.depth, mean(mean(radar.ages, 3), 2) + ...
    std(mean(radar.ages, 3), [], 2), 'm--')
plot(radar.depth, mean(mean(radar.ages, 3), 2) - ...
    std(mean(radar.ages, 3), [], 2), 'm--')
plot(cores.(name).depth, mean(cores.(name).ages,2), 'b', 'LineWidth', 2)
plot(cores.(name).depth, mean(cores.(name).ages,2) + ...
    std(cores.(name).ages, [], 2), 'b--')
plot(cores.(name).depth, mean(cores.(name).ages,2) - ...
    std(cores.(name).ages, [], 2), 'b--')



%%

for n = 1:20
    
    figure
    hold on
    i = randi(size(radar.ages,2));
    plot(radar.depth, ages_man(:,i), 'k', 'LineWidth', 2)
    plot(radar.depth, mean(radar.ages(:,i,:), 3), 'm', 'LineWidth', 2)
    plot(radar.depth, mean(radar.ages(:,i,:), 3) + ...
        std(squeeze(radar.ages(:,i,:)), [], 2), 'm--')
    plot(radar.depth, mean(radar.ages(:,i,:), 3) - ...
        std(squeeze(radar.ages(:,i,:)), [], 2), 'm--')
    plot(cores.(name).depth, mean(cores.(name).ages,2), 'b', 'LineWidth', 2)
    plot(cores.(name).depth, mean(cores.(name).ages,2) + ...
        std(cores.(name).ages, [], 2), 'b--')
    plot(cores.(name).depth, mean(cores.(name).ages,2) - ...
        std(cores.(name).ages, [], 2), 'b--')
    hold off
    
    figure
    hold on
    plot(radar.SMB_yr{i}, mean(radar.SMB{i}, 2), 'm', 'LineWidth', 2)
    plot(radar.SMB_yr{i}, mean(radar.SMB{i}, 2)+std(radar.SMB{i},[],2), 'm--')
    plot(radar.SMB_yr{i}, mean(radar.SMB{i}, 2)-std(radar.SMB{i},[],2), 'm--')
    plot(cores.(name).SMB_yr, mean(cores.(name).SMB, 2), 'b', 'LineWidth', 2)
    plot(cores.(name).SMB_yr, mean(cores.(name).SMB, 2) + ...
        std(cores.(name).SMB, [], 2), 'b--')
    plot(cores.(name).SMB_yr, mean(cores.(name).SMB, 2) - ...
        std(cores.(name).SMB, [], 2), 'b--')
    hold off
    
end
