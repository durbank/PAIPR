%  Function File: [p, F, df1, df2] = wanova (x1, x2)
%
%  Perform a Welch's alternative to one-way analysis of variance
%  (ANOVA). The goal is to test whether the population means of data
%  taken from k different groups are all equal. This test does not
%  require the condition of homogeneity of variances be satisfied.
%  For post-tests, it is recommended to run the function 'multicmp'.
%
%  Data should be given in a single vector y with groups specified by
%  a corresponding vector of group labels g (e.g., numbers from 1 to
%  k). This is the general form which does not impose any restriction
%  on the number of data in each group or the group labels.
%
%  Under the null of constant means, the Welch's test statistic F
%  follows an F distribution with df1 and df2 degrees of freedom.
%
%  The p-value (1 minus the CDF of this distribution at F) is
%  returned in the p output variable.
%
%  Bibliography:
%  [1] Welch (1951) On the Comparison of Several Mean Values: An
%       Alternative Approach. Biometrika. 38(3/4):330-336
%  [2] Tomarken and Serlin (1986) Comparison of ANOVA alternatives
%       under variance heterogeneity and specific noncentrality
%       structures. Psychological Bulletin. 99(1):90-99.
%
%  The syntax in this function code is compatible with most recent
%  versions of Octave and Matlab.
%
%  wanova v1.0 (last updated: 19/08/2013)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/

function [p, F, df1, df2] = wanova (varargin)

  if length(varargin)<2
    error('Invalid number of input arguments');
  end

  if nargout>4
    error('Invalid number of output arguments');
  end

%   if size(y,1)~=numel(y)
%     error('The first input argument must be a vector');
%   end
% 
%   if size(g,1)~=numel(g)
%     error('The second input argument must be a vector');
%   end

  % Determine the number of groups
  k = length(varargin);

  % Obtain the size, mean and variance for each sample
  n = zeros(k,1);
  mu = zeros(k,1);
  v = zeros(k,1);
  for i=1:k
    n(i,1) = numel(varargin{i});
    mu(i,1) = median(varargin{i});
%     mu(i,1) = mean(varargin{i});
    v(i,1) = var(varargin{i});
  end
  
  se = v./n;
  
  t = (mu - mu')./sqrt(se+se');
  nu = (se+se').^2./(se.^2./(n-1) + (se').^2./(n'-1));
  
  p = tcdf(t, nu);

  
  
  
  
  

%   % Take the reciprocal of the SEM
%   w = 1./se;
% 
%   % Calculate the origin
%   ori = sum(w.*mu)./sum(w);
% 
%   % Calculate Welch's test F ratio
%   F = (k-1)^-1*sum(w.*(mu-ori).^2)/...
%       (1+((2*(k-2)/(k^2-1))*sum((1-w/sum(w)).^2.*(n-1).^-1)));
% 
%   % Calculate the degrees of freedom
%   df1 = k-1;
%   df2 = (3/(k^2-1)*sum((1-w/sum(w)).^2.*(n-1).^-1))^-1;
% 
%   % Calculate the p-value
%   p = 1-fcdf(F,df1,df2);

end