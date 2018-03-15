% Function for stacking radar records by an arbitrary window size

function [mdata_stack] = radar_stack(mdata, window_length)

window_sz = round(window_length/mean(diff(mdata.dist)));

stack_idx = 1:window_sz:size(mdata.data_Z, 2);
E_stack = zeros(1, length(stack_idx)-1);
N_stack = zeros(1, length(stack_idx)-1);
data_stack = zeros(size(mdata.data_Z, 1), length(stack_idx)-1);
for i = 1:length(stack_idx)-1
    E_stack(i) = mean(mdata.Easting(stack_idx(i):stack_idx(i+1)));
    N_stack(i) = mean(mdata.Northing(stack_idx(i):stack_idx(i+1)));
    data_stack(:,i) = mean(mdata.data_Z(:,stack_idx(i):stack_idx(i+1)), 2);
end

dist_stack = 0:window_length:window_length*(length(stack_idx)-2);

mdata_stack = struct('Easting', N_stack, 'Northing', N_stack, 'dist', ...
    dist_stack, 'collect_date', mdata.collect_date, ...
    'depth', mdata.depth_interp, 'data_stack', data_stack);
end
    