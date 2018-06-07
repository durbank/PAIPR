% Function for stacking radar records by an arbitrary window size

function [mdata_stack] = radar_stack(mdata, window_length)

window_sz = round(window_length/mean(diff(mdata.dist)));

if size(mdata.depth, 2) > 1
    stack_idx = 1:window_sz:size(mdata.data_out, 2);
    E_stack = zeros(1, length(stack_idx)-1);
    N_stack = zeros(1, length(stack_idx)-1);
    data_stack = zeros(size(mdata.data_out, 1), length(stack_idx)-1);
    rho_stack = zeros(size(mdata.rho_coeff, 1), length(stack_idx)-1);
    var_stack = zeros(size(mdata.rho_var, 1), length(stack_idx)-1);
    depth_stack = zeros(size(mdata.data_out, 1), length(stack_idx)-1);
    
    for i = 1:length(stack_idx)-1
        E_stack(i) = mean(mdata.Easting(stack_idx(i):stack_idx(i+1)));
        N_stack(i) = mean(mdata.Northing(stack_idx(i):stack_idx(i+1)));
        data_stack(:,i) = mean(mdata.data_out(:,stack_idx(i):stack_idx(i+1)), 2);
        rho_stack(:,i) = mean(mdata.rho_coeff(:,stack_idx(i):stack_idx(i+1)), 2);
        var_stack(:,i) = mean(mdata.rho_var(:,stack_idx(i):stack_idx(i+1)), 2);
        depth_stack(:,i) = mean(mdata.depth(:,stack_idx(i):stack_idx(i+1)), 2);
    end
else
    stack_idx = 1:window_sz:size(mdata.data_out, 2);
    E_stack = zeros(1, length(stack_idx)-1);
    N_stack = zeros(1, length(stack_idx)-1);
    data_stack = zeros(size(mdata.data_out, 1), length(stack_idx)-1);
    for i = 1:length(stack_idx)-1
        E_stack(i) = mean(mdata.Easting(stack_idx(i):stack_idx(i+1)));
        N_stack(i) = mean(mdata.Northing(stack_idx(i):stack_idx(i+1)));
        data_stack(:,i) = mean(mdata.data_out(:,stack_idx(i):stack_idx(i+1)), 2);
    end
    rho_stack = repmat(mdata.rho_coeff, 1, size(data_stack, 2));
    var_stack = repmat(mdata.rho_var, 1, size(data_stack, 2));
    depth_stack = repmat(mdata.depth, 1, size(data_stack, 2));
end

dist_stack = 0:window_length:window_length*(length(stack_idx)-2);

if isfield(mdata, 'arr_layers')
    
    man_layers = mdata.arr_layers(:,stack_idx(1:end-1));
    mdata_stack = struct('collect_date', mdata.collect_date, 'Easting', E_stack,...
        'Northing', N_stack, 'dist', dist_stack,  'depth', depth_stack, ...
        'data_stack', data_stack, 'rho_coeff', rho_stack, 'rho_var', var_stack,...
        'man_layers', man_layers);
else
    
    mdata_stack = struct('collect_date', mdata.collect_date, 'Easting', E_stack,...
        'Northing', N_stack, 'dist', dist_stack,  'depth', depth_stack, ...
        'data_stack', data_stack, 'rho_coeff', rho_stack, 'rho_var', var_stack);
end
end
