figure
imagesc(peaks)
hold on
for i = 1:length(layers)
    [r,c] = ind2sub(size(peaks), layers{i});
    plot(c, r, '.', 'MarkerSize', 15)
end
hold off

figure
imagesc(radar.data_smooth, [-2 2])
hold on
for i = 1:length(layers_idx)
    [r,c] = ind2sub(size(peaks), layers_idx{i});
    plot(c, r, '.', 'MarkerSize', 15)
end
hold off