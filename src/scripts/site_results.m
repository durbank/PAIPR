% Script compare results from auto-IceBridge to cores and manual IceBridge

% Directories to data of interest based on computer (eventually will be
% replaced with GUI for data directory selection)
PC_true = ispc;
switch PC_true
    case true
%         computer = 'work';
        computer = input('Current PC: ');
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
name = 'SEAT10_4';

% Load relevant radar data (previously generated using the above section)
radar_OIB = load(fullfile(data_path, 'IceBridge/manual_layers', name, ...
    strcat('layers_', name, '_opt.mat')));

% Load manually traced layers for current SEAT core site (generated using
% 'draw_manual.m')
tmp = load(fullfile(data_path, 'IceBridge/manual_layers', name, ...
    strcat('manual_', name, '.mat')));
man_layers = tmp.man_all;

% Remove manual count layer data deeper than the deepest point in processed
% radargram data
keep_idx = cellfun(@(x) round(x(:,2))<=size(radar_OIB.data_smooth,1), ...
    man_layers, 'UniformOutput', false);
man_layers = cellfun(@(x,y) x(y,:), man_layers, keep_idx, 'UniformOutput', false);

%% Age-depth scales from manual count data

% Preallocate matrix array for manual layer group number indices
man_Gnum = zeros(size(radar_OIB.data_smooth));

% Place manual layer group numbers within a matrix array of the same size
% as radargram data
for i = 1:length(man_layers)
    
    layer_idx = sub2ind(size(man_Gnum), round(man_layers{i}(:,2)), ...
        man_layers{i}(:,1));
    man_Gnum(layer_idx) = i;
end

% Logical array of presence of manual layer member
man_log = logical(man_Gnum);

% Preallocate matrix of manual layer ages
ages_man = zeros(size(radar_OIB.data_smooth));

% Define surface age as the collection date of the radar data
age_top = radar_OIB.collect_date;
yr_pick1 = ceil(radar_OIB.collect_date - 1);

% Interpolate manual layer age-depth scales for each trace
for i = 1:size(ages_man,2)
    depths_i = [0; radar_OIB.depth(man_log(:,i))];
    yrs_i = ([age_top yr_pick1:-1:yr_pick1-length(depths_i)+2])';
    ages_man(:,i) = interp1(depths_i, yrs_i, radar_OIB.depth, 'linear', 'extrap');
end

%%

text_name = strrep(name, '_', '-');

% Find OIB data within 5 km of the ith core and the nearest trace to
% the core
dist_OIB = pdist2([cores.(name).Easting  cores.(name).Northing], ...
    [radar_OIB.Easting', radar_OIB.Northing']);
trace_idx = 1:length(radar_OIB.SMB);
OIB_idx = trace_idx(dist_OIB<=5000);
[~, OIB_near] = min(dist_OIB);

% Plot of age-depth estimates from OIB, SEAT core site, manual counts,
% and core data (only nearest trace)
f1 = figure('Position', [200 200 700 700]);
hold on
for n = 1:Ndraw
    h0 = plot(radar_OIB.depth, radar_OIB.ages(:,OIB_near,n), 'm', 'LineWidth', 0.5);
    h0.Color(4) = 0.05;
end
for n = 1:Ndraw
    h0 = plot(cores.(name).depth, cores.(name).ages(:,n), 'b', 'LineWidth', 0.5);
    h0.Color(4) = 0.05;
end
% h1 = plot(radar_i.depth, median(radar_i.ages(:,SEATi_near,:), 3),...
%     'r', 'LineWidth', 2);
% %     plot(radar_i.depth, median(radar_i.ages(:,SEATi_near,:), 3)...
% %         + 2*std(squeeze(radar_i.ages(:,SEATi_near,:)), [], 2), 'r--')
% %     plot(radar_i.depth, median(radar_i.ages(:,SEATi_near,:), 3)...
% %         - 2*std(squeeze(radar_i.ages(:,SEATi_near,:)), [], 2), 'r--')
h2 = plot(radar_OIB.depth, mean(radar_OIB.ages(:,OIB_near,:), 3),...
    'm', 'LineWidth', 2);
%     plot(radar_OIB.depth, median(radar_OIB.ages(:,OIB_near,:), 3)...
%         + 2*std(squeeze(radar_OIB.ages(:,OIB_near,:)), [], 2), 'm--')
%     plot(radar_OIB.depth, median(radar_OIB.ages(:,OIB_near,:), 3)...
%         - 2*std(squeeze(radar_OIB.ages(:,OIB_near,:)), [], 2), 'm--')
h3 = plot(cores.(name).depth, mean(cores.(name).ages, 2), 'b', 'LineWidth', 2);
%     plot(cores.(name).depth, mean(cores.(name).ages, 2) + ...
%         2*std(cores.(name).ages, [], 2), 'b--')
%     plot(cores.(name).depth, mean(cores.(name).ages, 2) - ...
%         2*std(cores.(name).ages, [], 2), 'b--')
h4 = plot(radar_OIB.depth, ages_man(:,OIB_near), 'k', 'LineWidth', 2);
legend([h2 h3 h4], 'Radar (auto)', 'Firn core', 'radar (manual)');
title(strcat(text_name, ' age-depth scale (nearest)'))
ylabel('Calendar years')
xlabel('Depth (m)')
hold off

f2 = figure('Position', [200 200 700 700]);
hold on
h1 = plot(radar_OIB.depth, mean(ages_man,2), 'k', 'LineWidth', 2);
plot(radar_OIB.depth, mean(ages_man,2) + std(ages_man, [], 2), ...
    'k--', 'LineWidth', 0.5)
plot(radar_OIB.depth, mean(ages_man,2) - std(ages_man, [], 2), ...
    'k--', 'LineWidth', 0.5)
h2 = plot(radar_OIB.depth, mean(mean(radar_OIB.ages, 3), 2), ...
    'm', 'LineWidth', 2);
plot(radar_OIB.depth, mean(mean(radar_OIB.ages, 3), 2) + ...
    std(mean(radar_OIB.ages, 3), [], 2), 'm--', 'LineWidth', 0.5)
plot(radar_OIB.depth, mean(mean(radar_OIB.ages, 3), 2) - ...
    std(mean(radar_OIB.ages, 3), [], 2), 'm--', 'LineWidth', 0.5)
h3 = plot(cores.(name).depth, mean(cores.(name).ages,2), ...
    'b', 'LineWidth', 2);
plot(cores.(name).depth, mean(cores.(name).ages,2) + ...
    std(cores.(name).ages, [], 2), 'b--', 'LineWidth', 0.5)
plot(cores.(name).depth, mean(cores.(name).ages,2) - ...
    std(cores.(name).ages, [], 2), 'b--', 'LineWidth', 0.5)
legend([h2 h3 h4], 'Radar (auto)', 'Firn core', 'radar (manual)');
title(strcat(text_name, ' age-depth scale (all)'))
ylabel('Calendar years')
xlabel('Depth (m)')
hold off

f3 = figure('Position', [200 200 1000 700]);
hold on
title(strcat(text_name, ' annual SMB (nearest)'))
for n = 1:Ndraw
    h0 = plot(radar_OIB.SMB_yr{OIB_near}, radar_OIB.SMB{OIB_near}(:,n), 'm', 'LineWidth', 0.5);
    h0.Color(4) = 0.05;
end
%         for n = 1:Ndraw
%             h0 = plot(radar_i.SMB_yr{SEATi_near}, radar_i.SMB{SEATi_near}(:,n), 'r', 'LineWidth', 0.5);
%             h0.Color(4) = 0.05;
%         end
for n = 1:Ndraw
    h0 = plot(cores.(name).SMB_yr, cores.(name).SMB(:,n), 'b', 'LineWidth', 0.5);
    h0.Color(4) = 0.05;
end
%     h1 = plot(radar_i.SMB_yr{SEATi_near}, median(radar_i.SMB{SEATi_near}, 2),...
%         'r', 'LineWidth', 2);
% %     plot(radar_i.SMB_yr{SEATi_near}, median(radar_i.SMB{SEATi_near}, 2) + ...
% %         2*std(radar_i.SMB{SEATi_near}, [], 2), 'r--');
% %     plot(radar_i.SMB_yr{SEATi_near}, median(radar_i.SMB{SEATi_near}, 2) - ...
% %         2*std(radar_i.SMB{SEATi_near}, [], 2), 'r--');
h2 = plot(radar_OIB.SMB_yr{OIB_near}, median(radar_OIB.SMB{OIB_near}, 2), ...
    'm', 'LineWidth', 2);
h3 = plot(cores.(name).SMB_yr, median(cores.(name).SMB, 2), ...
    'b', 'LineWidth', 2);
%     h4 = plot(man_yr, man_SMB(OIB_near), 'k', 'LineWidth', 2);

%     xlim([min([min(radar_i.SMB_yr{SEATi_near}) min(radar_OIB.SMB_yr{OIB_near}) ...
%         min(cores.(name).SMB_yr)]) ...
%         max([max(radar_i.SMB_yr{SEATi_near}) max(radar_OIB.SMB_yr{OIB_near}) ...
%         max(cores.(name).SMB_yr)])])
legend([h2 h3], 'OIB traces', 'Firn core')
xlabel('Calendar years')
ylabel('Annual SMB (mm w.e.)')
hold off


min_length = min(cellfun(@length, radar_OIB.SMB_yr));
SMB_OIB = cell2mat(cellfun(@(x) mean(x(1:min_length,:), 2), radar_OIB.SMB, ...
    'UniformOutput', false));
OIB_yr = radar_OIB.SMB_yr{1}(1:min_length);
f4 = figure('Position', [200 200 1000 700]);
hold on
title(strcat(text_name, ' annual SMB (all)'))
h1 = plot(OIB_yr, mean(SMB_OIB, 2), 'm', 'LineWidth', 2);
plot(OIB_yr, mean(SMB_OIB,2) + std(SMB_OIB,[],2), 'm--', 'LineWidth', 0.5)
plot(OIB_yr, mean(SMB_OIB,2) - std(SMB_OIB,[],2), 'm--', 'LineWidth', 0.5)
h2 = plot(cores.(name).SMB_yr, mean(cores.(name).SMB, 2), ...
    'b', 'LineWidth', 2);
plot(cores.(name).SMB_yr, mean(cores.(name).SMB, 2) + ...
    std(cores.(name).SMB,[],2), 'b--', 'LineWidth', 0.5)
plot(cores.(name).SMB_yr, mean(cores.(name).SMB, 2) - ...
    std(cores.(name).SMB,[],2), 'b--', 'LineWidth', 0.5)
legend([h1 h2], 'OIB traces', 'Firn core')
xlabel('Calendar years')
ylabel('Annual SMB (mm w.e.)')
hold off