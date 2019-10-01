% Function to optimize the logistic regression parameters using age-depth
% profiles from manually and autonomously traced layers

function [r_param, k_param, SSE] = opt_param(manual_vector, DB_vector, depth_vector)

% Contruct manual age-depth profile
yr_idx = find(manual_vector);
man_depth = [0; depth_vector(yr_idx)];
man_age = 0:length(yr_idx);


PAIPR_depth = [0; depth_vector(logical(DB_vector))];
age_interp = interp1(man_depth, man_age, PAIPR_depth);


keep_idx = ~isnan(age_interp);
PAIPR_depth = PAIPR_depth(keep_idx);
age_interp = age_interp(keep_idx);


null_val = -999;
DB_vals = [null_val; DB_vector(logical(DB_vector))];
DB_vals = DB_vals(keep_idx);
DB_vals(DB_vals == null_val) = [];


options = optimset('TolX',1e-2, 'MaxFunEvals', 25);
% options = optimset('PlotFcns','optimplotfval','TolX',1e-2, 'MaxFunEvals', 25);
x0 = [-3 3];
[xmin, fval] = fminsearch(@(x) ...
    loss_fun(PAIPR_depth,age_interp,DB_vals,x(1),x(2), 10000), x0, options);

r_param = xmin(1);
k_param = xmin(2);
SSE = fval;







end