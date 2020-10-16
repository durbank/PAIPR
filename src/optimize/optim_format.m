% Function to format interim radar data and manually-traced layers so they
% accurately overlap (stops biasing of results due to mismatched inputs)

function [radar] = optim_format(radar, man_layers)

% Get dimensions of overlapping results between manual and radar
man_sz = [round(max(cellfun(@(x) max(x(:,2)), man_layers))) ...
    round(max(cellfun(@(x) max(x(:,1)), man_layers)))];
rad_sz = size(radar.data_smooth);
dims = [min([man_sz(1) rad_sz(1)]) min([man_sz(2) rad_sz(2)])];

% Clip results to overlapping sections
man_keep = cellfun(...
    @(x) round(x(:,1)) <= dims(2) & round(x(:,2)) <= dims(1), ...
    man_layers, 'UniformOutput', false);
man_layers = cellfun(@(x,y) x(y,:), man_layers, man_keep, ...
    'UniformOutput', false);

% Subset man_layers to those that extend the full echogram length
last_idx = find(cellfun(@length, man_layers) >= dims(2), ...
    1, 'last');
layers_full = man_layers(1:last_idx);

% Create logical matrix of manual layers
man_idx = cellfun(@(x) sub2ind(dims, round(x(:,2)),x(:,1)), ...
    layers_full, 'UniformOutput', false);
man_grid = zeros(dims);
for k = 1:length(man_idx)
    man_grid(man_idx{k}) = k;
end
radar.man_grid = man_grid;


% Clip remaining radar fields to match overlapping section
flds = fieldnames(radar);
for i=1:length(flds)
    
    fld = radar.(flds{i});
    dim_i = size(fld);
    row_idx = min([dims(1) dim_i(1)]);
    col_idx = min([dims(2) dim_i(2)]);
    radar.(flds{i}) = fld(1:row_idx,1:col_idx);
end

% Add clipped (but not subsetted) manual layers (with indices rather than 
% row/col) to radar structure
null_idx = cellfun(@isempty, man_layers);
man_idx = cellfun(@(x) sub2ind(dims, round(x(:,2)),x(:,1)), ...
    man_layers(~null_idx), 'UniformOutput', false);
radar.man_layers = man_idx;

end