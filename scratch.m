figure
imagesc(peaks_raw)
hold on
for i = 1:length(layers)
    [r,c] = ind2sub(size(peaks_raw), layers{i});
    plot(c, r, '.', 'MarkerSize', 15)
end
hold off

figure
imagesc(radar.data_smooth, [-2 2])
hold on
for i = 1:length(radar.layers)
    [r,c] = ind2sub(size(radar.data_smooth), radar.layers{i});
    plot(c, r, 'LineWidth', 3)
end
hold off