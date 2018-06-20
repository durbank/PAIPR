figure
imagesc(peaks_raw)
hold on
for k = 1:length(layers)
    [r,c] = ind2sub(size(peaks_raw), layers{k});
    plot(c, r, '.', 'MarkerSize', 15)
end
hold off

figure
imagesc(radar.data_smooth, [-2 2])
hold on
for k = 1:length(radar.layers)
    [r,c] = ind2sub(size(radar.data_smooth), radar.layers{k});
    plot(c, r, 'LineWidth', 3)
end
hold off


%% 
% Trace idx to investigate
i = randi(size(radar.data_smooth, 2));

% Find the nearest cores to the radar data (for comparison plots)
[~, cores_near_idx] = sort(pdist2([radar.Easting(i) radar.Northing(i)], ...
    [cores.Easting' cores.Northing'], 'Euclidean'));
core_near1 = cores.(cores.name{cores_near_idx(1)});
core_near2 = cores.(cores.name{cores_near_idx(2)});
core_near3 = cores.(cores.name{cores_near_idx(3)});

% Calculate the mean age-depth scale and std for radar trace i
age_mean = mean(squeeze(radar.ages(:,i,:)), 2);
age_std = std(squeeze(radar.ages(:,i,:)), [], 2);

% Plot full radargram
yr_idx = logical([diff(floor(age_mean)); 0]);
depth = radar.depth(yr_idx);
col = i*ones(length(depth),1);
% row = find(radar.likelihood(:,i)>0.75);
% col = i*ones(length(row),1);
figure('Position', [200 200 1500 800])
imagesc(radar.dist, radar.depth, radar.data_smooth, [-2 2])
colorbar
xlabel('Distance along profile (m)')
ylabel('Depth (m)')
hold on
plot(radar.dist(col), depth, 'r.', 'MarkerSize', 25)
xlim([0 radar.dist(end)])
ylim([0 radar.depth(end)])
set(gca, 'Ydir', 'reverse', 'FontSize', 18)
hold off

% Age-depth scale comparison between radar trace and nearest cores
figure
hold on
h1 = plot(core_near1.depth, mean(core_near1.ages, 2), 'b', 'LineWidth', 2);
plot(core_near1.depth, mean(core_near1.ages, 2) + 2*std(core_near1.ages, [], 2), 'b--')
plot(core_near1.depth, mean(core_near1.ages, 2) - 2*std(core_near1.ages, [], 2), 'b--')
h2 = plot(core_near2.depth, core_near2.age, 'c', 'LineWidth', 2);
h3 = plot(core_near3.depth, core_near3.age, 'c--', 'LineWidth', 1);
h4 = plot(radar.depth, age_mean, 'r', 'LineWidth', 2);
plot(radar.depth, age_mean + 2*age_std, 'r--', 'LineWidth', 0.5)
plot(radar.depth, age_mean - 2*age_std, 'r--', 'LineWidth', 0.5)
ylabel('Calendar Year')
xlabel('Depth (m)')
legend([h1 h2 h3 h4], 'Nearest core age (manual)', '2nd nearest core', ...
    '3rd nearest core', 'Radar age (automated)', 'Location', 'ne')
set(gca, 'FontSize', 10)
hold off
