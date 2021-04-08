% Script to compare logistic parameters between different echograms and
% calculated bulk parameter statistics

% Set environment
DATA_DIR = ['/media/durbank/WARP/Research/Antarctica/Data/'...
    'IceBridge/optimization/v0.4.0/'];

% Get list of param files
files = dir(fullfile(DATA_DIR, 'params', '*.mat'));
str_rm = 'params_flight_';

% Preallocation
params = struct();
xi = -16:0.05:0; % Limits distribution to reasonable values
ks_dens = zeros(length(xi),length(files));

for i=1:length(files)
    
    % Load data and create more useful variable names
    data = load(fullfile(files(i).folder, files(i).name));
    f_parts = split(extractAfter(files(i).name, str_rm), '_');
    tmp = strjoin(f_parts(2:end), '_');
    f1 = split(tmp, '.');
    f_name = strrep(strcat(f1{1}, '_', f_parts{1}), '-', '_');
    
    % Assign values to output structure
    params.(f_name).r_all = data.r_params;
    params.(f_name).SSE_all = data.SSE;
    
    % Calculate weights and r statistics for each file
    weights = (1./data.SSE) ./ sum(1./data.SSE);
    params.(f_name).r_mu = sum(weights .* data.r_params);
    params.(f_name).r_std = std(data.r_params, weights);
    
    % Get estimated SSE-weighted distribution for r parameter associated
    % with each file
    ks_dens(:,i) = ksdensity(data.r_params, xi, 'Weights', weights);
end

% Get list of file fieldnames
p_flds = fieldnames(params);

% Figure comparing the various distributions in r values for each file
% (color-matched for overlap between 2011 and 2016)
figure('Position', [10 10 900 500])
hold on
plot(xi, ks_dens(:,1), 'Color', [0.9290, 0.6940, 0.1250])
plot(xi, ks_dens(:,2), 'r')
plot(xi, ks_dens(:,3), 'cyan')
plot(xi, ks_dens(:,4), 'black')
plot(xi, ks_dens(:,5), 'b')
plot(xi, ks_dens(:,6), 'm')
plot(xi, ks_dens(:,7), 'r--')
plot(xi, ks_dens(:,8), 'b--')
plot(xi, ks_dens(:,9), 'm--')
legend(p_flds, 'Interpreter', 'none', 'Location', 'northwest')
hold off

% Calculate the bulk SSE-weighted r-value and r st dev.
tmp = cellfun(@(x) x.r_all, struct2cell(params), 'UniformOutput', false);
r_all = [tmp{:}];
tmp = cellfun(@(x) x.SSE_all, struct2cell(params), 'UniformOutput', false);
SSE_all = [tmp{:}];
weights = (1./SSE_all) ./ sum(1./SSE_all);

% Print bulk r statistics to command window
sprintf('Mean SSE-weighted r value for all training sets: %.3g', ...
    sum(weights .* r_all))
sprintf('SSE-weighted st. dev. of r value for all training sets: %.3g',...
    std(r_all, weights))

    