% Script to generate a mean density variance with depth model based on the
% spatially weighted density variances of nearby cores

function [rho_std] = rho_variance(cores, core_idx, SWM)


core_names = cores.name(core_idx);
rho_std = cell(1, length(core_names));

for i = 1:length(core_names)
    core_i = cores.(core_names{i});
    rho_std{i} = sqrt(movvar(core_i.rho, round(length(rho_res)/10)));
end



end