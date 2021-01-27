% Function to calculate layer likelihood scores and age-depth profile
% distributions (based on likelihood scores)

function [radar, Ndraw_new] = radar_age(radar, r_bar, r_sigma, k, Ndraw)
%% Assign layer likelihood scores and estimate age-depth scales


% N_sqrt = floor(sqrt(Ndraw));
% N_sqrt2 = N_sqrt;

% % Define normal distribution of r-parameter for MC simulations
% r_draws = normrnd(r_bar, r_sigma, 1, N_sqrt); %This may need to move within the loop to save memory

% Preallocate arrays for layer likelihoods and ages
% err_out = [];
ages = NaN([size(radar.data_smooth) Ndraw]);
% ages = zeros([size(radar.data_smooth) N_sqrt N_sqrt2]);
% like_all = zeros([size(radar.data_smooth) N_sqrt]);


for i = 1:size(radar.DB, 2)
    % Define surface age and the year associated with the first pick 
    % of the algorithm
    yr_vec = datevec(radar.collect_time(i));
    yr_pick1 = yr_vec(1);
    age_top = yr_vec(1) + (30*yr_vec(2)+yr_vec(3))/365;
    
    % Get layer prom-distance values and depths for layers in ith trace
    peaks_i = radar.DB(:,i);
    peaks_idx = peaks_i>0;
    peaks_i = peaks_i(peaks_idx);
    depths_i = radar.depth(peaks_idx);
    
    
%     likelihood = 1./(1+exp(r_draws.*peaks_i + k));
%     like_all(peaks_idx,i,:) = likelihood;
    
    for MC=1:Ndraw
        
        % Sample r (this might be better modeled as a beta distribution?) 
        % for iteration and calculate likelihoods
        r = normrnd(r_bar, r_sigma); %
        likelihood = 1./(1+exp(r.*peaks_i + k));
        
        % Generate logical vector for which layers should be years in this
        % MC iteration
        yr_idx = zeros(length(depths_i));
        for j=1:length(yr_idx)
            yr_idx(j) = rand(1) <= likelihood(j);
        end
        
        % Get vector of depths and integer ages for this MC iteration
        depths_MC = [0; depths_i(logical(yr_idx))];
        yrs_MC = ([age_top yr_pick1:-1:yr_pick1-length(depths_MC)+2])';
        
        try
            ages(:,i,MC) = interp1(depths_MC, yrs_MC, ...
                radar.depth, 'linear', 'extrap');
        catch
            sprintf(['Error in age interpolation for trace %u, '...
                'trial %u. Filling with mean ages.'], i, j)
%             err_out = [err_out MC];
        end
    end
%     if ~isempty(err_out)
%         ages(:,i,err_out) = repmat(sum(squeeze(ages(:,i,:)), 2)./...
%             sum(squeeze(ages(:,i,:))~=0, 2), 1, length(err_out));
%     end
%     err_out = [];
        
end
    
    
%         % Assign MC simulation annual layer presence based on layer likelihood
%         % values
%         yr_idx = zeros(length(depths_i), N_sqrt2);
%         for j = 1:length(depths_i)
%             R = rand(N_sqrt2, 1) <= likelihood(j,MC);
%             yr_idx(j,:) = R;
%         end
%         
%         for j = 1:N_sqrt2
%             
%             depths_j = [0; depths_i(logical(yr_idx(:,j)))];
%             yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
%             try
%                 ages(:,i,MC,j) = interp1(depths_j, yrs_j, radar.depth, ...
%                     'linear', 'extrap');
%             catch
%                 sprintf(['Error in age interpolation for trace %u, '...
%                     'trial %u. Filling with mean ages.'], i, j)
%                 err_out = [err_out j];
%             end
%         end
%         if ~isempty(err_out)
%             ages(:,i,MC,err_out) = repmat(sum(squeeze(ages(:,i,MC,:)), 2)./...
%                 sum(squeeze(ages(:,i,MC,:))~=0, 2), 1, length(err_out));
%         end
%         err_out = [];
%     end
% end
% 
% 
% ages = reshape(ages, [size(radar.data_smooth) N_sqrt*N_sqrt2]);
% 
% % Linearly interpolate traces with missing or problematic age-depth scales
% radar.ages = fillmissing(ages, 'linear', 2);
% 
% 
% radar.likelihood = like_all;
% Ndraw_new = N_sqrt*N_sqrt2;

% Fill in missing age results via interpolation with surrounding traces
radar.ages = fillmissing(ages, 'linear', 2);
% radar.ages = ages;

% Generate likelihood image based simply on r_bar
likelihood = zeros(size(radar.data_smooth));
likelihood(radar.DB > 0) = 1./(1+exp(r_bar.*radar.DB(radar.DB > 0) + k));
radar.likelihood = likelihood;

% For legacy support, assign Ndraw_new
Ndraw_new = Ndraw;














% for MC=1:length(N_sqrt)
%     
%     for i = 1:size(radar.DB, 2)
%         
%         % Define surface age and the year associated with the first pick 
%         % of the algorithm
%         yr_vec = datevec(radar.collect_time(i));
%         yr_pick1 = yr_vec(1);
%         age_top = yr_vec(1) + (30*yr_vec(2)+yr_vec(3))/365;
%         
%         % Get layer prom-distance values and depths for layers in ith trace
%         peaks_i = radar.DB(:,i);
%         peaks_idx = peaks_i>0;
%         peaks_i = peaks_i(peaks_idx);
%         depths_i = radar.depth(peaks_idx);
%         
%         %     likelihood = 1./(1+exp(r*peaks_i + k));
%         %     radar.likelihood(peaks_idx,i) = likelihood;
%         likelihood = 1./(1+exp(r_draws(MC)*peaks_i + k));
%         like_all(peaks_idx,i, MC) = likelihood;
%         
%         % Assign MC simulation annual layer presence based on layer likelihood
%         % values
%         yr_idx = zeros(length(depths_i), N_sqrt);
%         for j = 1:length(depths_i)
%             R = rand(N_sqrt, 1) <= likelihood(j);
%             yr_idx(j,:) = R;
%         end
%         
%         for j = 1:N_sqrt
%             
%             depths_j = [0; depths_i(logical(yr_idx(:,j)))];
%             yrs_j = ([age_top yr_pick1:-1:yr_pick1-length(depths_j)+2])';
%             try
%                 ages(:,i,j) = interp1(depths_j, yrs_j, radar.depth, ...
%                     'linear', 'extrap');
%             catch
%                 sprintf('Error in age interpolation for trace %u, trial %u. Filling with mean ages.', i, j)
%                 err_out = [err_out j];
%             end
%         end
%         if ~isempty(err_out)
%             ages(:,i,err_out) = repmat(sum(squeeze(ages(:,i,:)), 2)./...
%                 sum(squeeze(ages(:,i,:))~=0, 2), 1, length(err_out));
%         end
%         err_out = [];
%     end
% end

end