function [h] = EstimateBWDensCV(X)

  % Assumes gaussian kernel
  N      = length(X);
  Kstar  = @(u) normpdf(u,0,sqrt(2)) - 2*normpdf(u);
  K0     = normpdf(0);
  objfcn = @(h) sum(sum( Kstar( (repmat(X,1,N)-repmat(X',N,1))/h ) ))/((N^2)*h)...
                + 2*K0/(N*h);

  h = fminbnd(objfcn, eps(), 100);
end

