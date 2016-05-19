function [h] = EstimateBWRegCV(X,Y)

  function [err] = ObjectiveFunction(N,X,Y,h)
    r = nan(N,1);
    for n = 1:N
      use  = [1:n-1 n+1:N];
      r(n) = NadayaraWatson(Y(use),X(use),X(n),h);
    end
    err = sum((Y - r).^2)/N;
  end

  N = length(X);
  objfcn = @(h) ObjectiveFunction(N,X,Y,h);
  h = fminbnd(objfcn, eps(), 100);
end

