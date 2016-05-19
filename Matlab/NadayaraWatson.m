function [r] = NadayaraWatson(Y,X,x,h)
  N  = length(X); % Number of observations
  M  = length(x); % Number of evaluation points

  % Use Gaussian Kernel
  K = @(x) normpdf(x);

  Ks = K( (repmat(X,1,M)-repmat(x',N,1))/h );
  r  = [(Y'*Ks)./sum(Ks)]';
end

