%% OVERVIEW - mvnrnd.m
%
% Function to draw from an MVn random variable
function [R] = mvnrnd(MU,SIGMA,cases)

if nargin < 3
  cases=1;
end
if size(MU,2) > 1
  MU = MU'
end

nel = length(MU);
L   = chol(SIGMA)';

R = repmat(MU,1,cases) + L*randn(nel,cases);

end
