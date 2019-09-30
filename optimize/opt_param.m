% Function to optimize the logistic regression parameters using age-depth
% profiles from manually and autonomously traced layers

function [r_param, k_param] = opt_param(manual_vector, DB_vector, depth_vector)

% Contruct manual age-depth profile
yr_idx = find(manual_vector);
man_depths = depth_vector(yr_idx);
man_ages = 1:length(yr_idx);


PPR_depths = depth_vector(logical(DB_vector));

r_param = 0;

end