% Script to generate a mean density variance with depth model based on the
% spatially weighted density variances of nearby cores

function [var_param] = rho_variance(cores, core_idx, SWM)

core_names = cores.name(core_idx);
var_coeff = zeros(3, length(core_names));
for i = 1:length(core_names)
    core_i = cores.(core_names{i});
    rho_var = movvar(core_i.rho, round(length(core_i.rho)/5));
    
    fun = @(a, x) (mean(rho_var(1:25))-mean(rho_var(end-25:end)))./(a.^x) + ...
        mean(rho_var(end-25:end));
    
    var_coeff(1,i) = mean(rho_var(1:25));
    var_coeff(2,i) = mean(rho_var(end-25:end));
    var_coeff(3,i) = nlinfit(core_i.depth, rho_var, fun, 1.5);
    
    %     rho_top = mean(core_i.rho(1:10));
    %     g = @(r, x) (rho_ice-rho_top)*(x)./(x+r) + rho_top;
    %     r(i) = nlinfit(core_i.depth, core_i.rho, g, 20);
    
end

var_param = sum(SWM.*var_coeff, 2);

end