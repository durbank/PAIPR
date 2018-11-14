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
SMB = [mSMB2 mSMB];
b = [beta2 beta];

figure
imagesc(E, radar.depth, data, [-2 2])

figure
yyaxis left
plot(E, mean(SMB), 'b')
yyaxis right
plot(E, b, 'r')

figure
imagesc(radar2.Easting, radar2.depth, radar2.data_smooth, [-2 2])

figure
yyaxis left
plot(radar2.Easting, mean(mSMB2), 'b')
yyaxis right
plot(radar2.Easting, beta2, 'r')

figure
hold on
plot(radar.Easting, beta)
plot(radar2.Easting, beta2)
hold off

