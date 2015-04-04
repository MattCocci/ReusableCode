% OVERVIEW - mvnconddist.m
%
% Suppose X and Y that are jointly MVN such that
% 
% | Y | ~ MVN( |mu_Y|, |S_YY S_YX| )
% | X |      ( |mu_X|  |S_XY S_XX| )
%
% Then having observed X, this function will return the mean and the
% variance for the conditional distribution Y|X, also know as the
% regression estimate of Y given X
function [varargout] = mvnconddist(X, mu_Y, mu_X, S_YY, S_XX, S_YX)

  YbarX.mu    = mu_Y + S_YX*S_XX\(X-mu_X);
  YbarX.Sigma = S_YY - S_YX*S_XX\(S_YX');
  
  if nargout > 1
    varargout{1} = YbarX.mu;
    varargout{2} = YbarX.Sigma; 
  else
    varargout{1} = YbarX; 
  end
end
