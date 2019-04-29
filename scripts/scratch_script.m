

%% Testing accum_distGamma.m

addpath('../')
radar = load('/media/durbank/WARP/Research/Antarctica/Data/IceBridge/SEAT10_4to10_6/2011_SNO/SMB_results/radar_out6.mat');

data_tables = cellfun(@accum_distGamma, radar.SMB, radar.SMB_yr, ...
    'UniformOutput', false);

min_length = min(cellfun(@(x) length(x.Year), data_tables));
accum_med = cell2mat(cellfun(@(x) x.Median(1:min_length), data_tables, ...
    'UniformOutput', false));


years = 2010:-1:(2010-min_length+1);

figure
boxplot(fliplr(accum_med'), 'Labels', fliplr(cellfun(@num2str, num2cell(years), ...
    'UniformOutput', false)))
hold on
line(length(years):-1:1, accum_med(:,randi(size(accum_med,2))), 'Color', 'red')
line(length(years):-1:1, accum_med(:,randi(size(accum_med,2))), 'Color', 'cyan')
line(length(years):-1:1, accum_med(:,randi(size(accum_med,2))), 'Color', 'green')
hold off

