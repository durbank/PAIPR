%% Why are there such stark differences in trends between adjacent
%  radargrams, even when the mean SMB doesn't have nearly as large as jump?



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

wild = '*.mat';
OIB_files = dir(fullfile(data_path, 'IceBridge/SNO_radar/',...
    '2011/SMB_results/', wild));

radar = load(fullfile(OIB_files(19).folder, OIB_files(19).name));
radar2 = load(fullfile(OIB_files(18).folder, OIB_files(18).name));

yrs = (2010:-1:2010-32+1)';
mSMB = cell2mat(cellfun(@(x) mean(x(1:32,:),2), radar.SMB, ...
    'UniformOutput', false));
mSMB2 = cell2mat(cellfun(@(x) mean(x(1:32,:),2), radar2.SMB, ...
    'UniformOutput', false));

beta = zeros(1,length(mSMB));
for i=1:length(mSMB)
    coeff = robustfit(yrs, mSMB(:,i));
    beta(i) = coeff(2);
end

beta2 = zeros(1,length(mSMB2));
for i=1:length(mSMB2)
    coeff = robustfit(yrs, mSMB2(:,i));
    beta2(i) = coeff(2);
end

E = [radar2.Easting radar.Easting];
data = [radar2.data_smooth radar.data_smooth];
% SMB = [mSMB2 mSMB];
% b = [beta2 beta];



layers = cell(1, max(radar.groups(:)));
for i=1:length(layers)
    [r,c] = find(radar.groups==i);
    layers{i} = [r c];
end
layers = layers(~cellfun(@isempty, layers));

layers2 = cell(1, max(radar2.groups(:)));
for i=1:length(layers2)
    [r,c] = find(radar2.groups==i);
    layers2{i} = [r c];
end
layers2 = layers2(~cellfun(@isempty, layers2));

figure
imagesc(E, radar.depth, data, [-2 2])
hold on
cellfun(@(x) plot(radar.Easting(x(:,2)), radar.depth(x(:,1)), ...
    'r', 'LineWidth', 1.5), layers);
cellfun(@(x) plot(radar2.Easting(x(:,2)), radar2.depth(x(:,1)), ...
    'm', 'LineWidth', 1.5), layers2);

figure
yyaxis left
hold on
plot(radar.Easting, mean(mSMB), 'b')
plot(radar2.Easting, mean(mSMB2), 'c')
yyaxis right
plot(radar.Easting, beta, 'r')
plot(radar2.Easting, beta2, 'm')
xlim([min(E) max(E)])
hold off

figure
hold on
num_plot = 100;
for n = 1:num_plot
    h0 = plot(yrs, mSMB(:,length(mSMB)-n), 'r', 'LineWidth', 0.5);
    h0.Color(4) = 0.05;
end
for n = 1:num_plot
    h0 = plot(yrs, mSMB2(:,length(mSMB2)-n), 'm', 'LineWidth', 0.5);
    h0.Color(4) = 0.05;
end
hold off


