% Script to import data from OIB_Ku, OIB_SNOW, and SEAT2010 for the
% SEAT2010-4 core site, and compare nearest trace points

data_dir = '/Volumes/WARP/Research/Antarctica/WAIS Variability/accum-radar/output/';

load(strcat(data_dir, 'SEAT10_4.mat'));
load(strcat(data_dir, 'SNOW10_4.mat'));
load(strcat(data_dir, 'KU10_4.mat'));

%%

SEAT_idx = 34;
E = SEAT10_4_SMB.Easting(SEAT_idx);
N = SEAT10_4_SMB.Northing(SEAT_idx);

%compute Euclidean distances:
distances = sqrt(sum(bsxfun(@minus, [KU10_4_SMB.Easting' KU10_4_SMB.Northing'],...
    [E N]).^2,2));

% Index of nearest trace
[~, KU_idx] = min(distances);

%compute Euclidean distances:
distances = sqrt(sum(bsxfun(@minus, [SNOW10_4_SMB.Easting' SNOW10_4_SMB.Northing'],...
    [E N]).^2,2));

% Index of nearest trace
[~, SNOW_idx] = min(distances);

figure
hold on
plot(SEAT10_4_SMB.radar_yr, SEAT10_4_SMB.radar_accum(:,SEAT_idx), 'b', 'LineWidth', 2)
plot(KU10_4_SMB.radar_yr, KU10_4_SMB.radar_accum(:,KU_idx), 'r', 'LineWidth', 2)
plot(SNOW10_4_SMB.radar_yr, SNOW10_4_SMB.radar_accum(:,SNOW_idx), 'm', 'LineWidth', 2)
plot(SEAT10_4_SMB.radar_yr, SEAT10_4_SMB.radar_accum(:,SEAT_idx) + ...
    SEAT10_4_SMB.radar_ERR(:,SEAT_idx), 'b--', 'LineWidth', 0.5)
plot(SEAT10_4_SMB.radar_yr, SEAT10_4_SMB.radar_accum(:,SEAT_idx) - ...
    SEAT10_4_SMB.radar_ERR(:,SEAT_idx), 'b--', 'LineWidth', 0.5)
plot(KU10_4_SMB.radar_yr, KU10_4_SMB.radar_accum(:,KU_idx) + ...
    KU10_4_SMB.radar_ERR(:,KU_idx), 'r--', 'LineWidth', 0.5)
plot(KU10_4_SMB.radar_yr, KU10_4_SMB.radar_accum(:,KU_idx) - ...
    KU10_4_SMB.radar_ERR(:,KU_idx), 'r--', 'LineWidth', 0.5)
plot(SNOW10_4_SMB.radar_yr, SNOW10_4_SMB.radar_accum(:,SNOW_idx) + ...
    SNOW10_4_SMB.radar_ERR(:,SNOW_idx), 'm--', 'LineWidth', 0.5)
plot(SNOW10_4_SMB.radar_yr, SNOW10_4_SMB.radar_accum(:,SNOW_idx) - ...
    SNOW10_4_SMB.radar_ERR(:,SNOW_idx), 'm--', 'LineWidth', 0.5)
legend('SEAT Ku', 'OIB Ku', 'OIB SNOW')
hold off
