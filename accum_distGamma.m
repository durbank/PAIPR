% Function to fit gamma distributions to annual accumulation results and
% output distribution parameters

% This takes Monte Carlo draws of annually-resolved time series of accumulation 
% results and returns fitted gamma distribution parameters for each year

% Data inputs should be an m x n matrix, where m is the number of years in
% the time series and n is the number of Monte Carlo simulations

function [accum_table] = accum_distGamma(accum_data, accum_yr)

gamma_params = zeros(size(accum_data,1), 6);

for i = 1:size(accum_data, 1)
    
    data_i = accum_data(i,:)';
    
    [phat, pci] = gamfit(data_i);
    % pd_gamma = fitdist(data_i, 'Gamma');
    
    gamma_params(i,1) = phat(1);
    gamma_params(i,2:3) = pci(:,1);
    gamma_params(i,4) = phat(2);
    gamma_params(i,5:6) = pci(:,2);
    
    
end

accum_table = table(accum_yr, median(accum_data, 2), gamma_params(:,1), ...
    gamma_params(:,4), gamma_params(:,2:3), gamma_params(:,5:6), 'VariableNames', ...
    {'Year', 'Median', 'gamma_shape', 'gamma_scale', 'Shape_CI95', 'Scale_CI95'});
    


end





