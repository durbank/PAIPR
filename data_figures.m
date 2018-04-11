% Script to produce maps of the location of imported Antarctic data and
% various other plots potentially useful in the future

PC_true = ispc;
switch PC_true
    case true
        data_path = 'E:\Research\Antarctica\WAIS Variability\';
        addon_path = 'E:\Research\Antarctica\WAIS Variability\Addons\';
    case false
        data_path = '/Volumes/WARP/Research/Antarctica/WAIS Variability/';
        addon_path = '/Users/Durbank/Documents/MATLAB/Add-Ons/';
end

% Addons needed for analysis
% Add Antarctic Mapping Toolbox (AMT) to path
addon_folder = strcat(addon_path, 'AntarcticMappingTools_v5.03/');
addpath(genpath(addon_folder))

[cores] = import_cores(strcat(data_path, ['SEAT_cores' filesep 'DGK_core_data.xlsx']));
[SEAT10_GPR] = concat_loc(strcat(data_path, ['SEAT_Traverses' filesep ...
    'SEAT2010Kuband' filesep 'RawSEAT2010' filesep]), 'layers_*');
[SEAT11_GPR] = concat_loc(strcat(data_path, ['SEAT_Traverses' filesep ...
    'SEAT2011Kuband' filesep 'RawSEAT2011' filesep]), 'layers_*');
labels = strrep(cores.name, '_', '-');
basins = shaperead(strcat(data_path, 'Data/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'));

Easting_lims = [min(cores.Easting)-0.10*(max(cores.Easting)-min(cores.Easting)) ...
    max(cores.Easting)+0.15*(max(cores.Easting)-min(cores.Easting))];
Northing_lims = [min(cores.Northing)-0.1*(max(cores.Northing)-min(cores.Northing)) ...
    max(cores.Northing)+0.1*(max(cores.Northing)-min(cores.Northing))];
elev = cryosat2_data(Easting_lims, Northing_lims);

figure('Position', [10 10 1400 800])
hold on
h0 = image(Easting_lims, Northing_lims, elev, 'CDataMapping', 'scaled');
colormap(gray)
mapshow(basins, 'FaceAlpha', 0, 'LineWidth', 3)
h1 = scatter(cores.Easting, cores.Northing, 100, 'b', 'filled');
h2 = plot(SEAT10_GPR.Easting(1), SEAT10_GPR.Northing(1), 'r', 'LineWidth', 2);     % Correctly display radar as line in legend
plot(SEAT10_GPR.Easting, SEAT10_GPR.Northing, 'r.', 'MarkerSize', 0.05)
plot(SEAT11_GPR.Easting, SEAT11_GPR.Northing, 'r.', 'MarkerSize', 0.05)
text(cores.Easting, cores.Northing, strcat('\leftarrow', labels), 'FontSize', 18, 'Interpreter', 'tex');
h3 = plot(cores.Easting(1), cores.Northing(1), 'k', 'LineWidth', 2);
c0 = colorbar;
c0.Label.String = 'Elevation (m asl)';
c0.Label.FontSize = 18;
graticuleps(-81:0.5:-77,-125:2:-105, 'c')
xlim(Easting_lims)
ylim(Northing_lims)
scalebarps
box on
mapzoomps('ne', 'insetsize', 0.30)
legend([h1 h2 h3], 'Firn cores', 'Radar transects', 'WAIS Divide', 'Location', 'northwest')
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 18)
hold off


% %%% Code from Clem for plotting annual layers in radar
% 
% % Load radar data and clean to standard format
% [radar] = radar_clean(['SEAT Traverses' filesep 'GPR_SEAT2010_4_to_5.mat']);
% 
% % Preallocate matrix
% layer_idx = zeros(max(max(radar.arr_layers)),length(radar.dist));
% 
% % For loop to find the indexed location of each picked horizon at each
% % traverse location
% for i = 1:max(max(radar.arr_layers))
% [row,col] = find(radar.arr_layers==i);
% layer_idx(i,col) = row;
% end
% clear row col i
% 
% % Change back 0 to NaN
% layer_idx(layer_idx==0) = NaN;
% 
% % Transpose matrix for ease of plotting
% layer_idx = layer_idx';
% 
% % Plot picked layers over the radar image
% figure 
% hold on
% imagesc(radar.smooth)
% plot(layer_idx)
% set(gca, 'Ydir', 'reverse')
% hold off