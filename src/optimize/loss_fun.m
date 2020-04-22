% Loss function (based on SSE for age-depth residuals between manually and
% autonomously traced layers) used to optimize logistic regression
% parameters

function [f] = loss_fun(depth, age_interp, L_values, r, k, Ndraw)

ages = zeros(length(depth), Ndraw);
likelihood = [1; 1./(1+exp(r*L_values + k))];

% Assign MC simulation annual layer presence based on layer likelihood
% values
yr_idx = zeros(length(depth), Ndraw);
for j = 1:length(depth)
    R = rand(Ndraw, 1) <= likelihood(j);
    yr_idx(j,:) = R;
end


% for j = 1:Ndraw
%     depths_j = depth(logical(yr_idx(:,j)));
%     yrs_j = (0:(length(depths_j)-1))';
%     if length(depths_j) < 2
%         ages(:,j) = interp1([0 depth(end)], [0 1], depth, ...
%             'linear', 'extrap');
%     else
%         ages(:,j) = interp1(depths_j, yrs_j, depth, 'linear', 'extrap');
%     end
% end

for j = 1:Ndraw
    depths_j = depth(logical(yr_idx(:,j)));
    yrs_j = (0:(length(depths_j)-1))';
    try
        ages(:,j) = interp1(depths_j, yrs_j, depth, 'linear', 'extrap');
    catch
        ages(:,j) = NaN;
    end
end

age_out = interp1(depth, nanmean(ages,2), depth);


f = sum((age_out - age_interp).^2);



end