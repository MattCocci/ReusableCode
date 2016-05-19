function [reg] = Regress(data, Yname, Xnames)

  N = length(data.(Yname));
  K = length(Xnames)+1;

  % Get the Y matrix
  Y = data.(Yname);

  % Build the X matrix
  X = nan(N,K);
  X(:,1) = 1;
  for k = 2:K
    X(:,k) = data.(Xnames{k-1});
  end

  % OLS Results
  reg.Y     = Y;
  reg.X     = X;
  reg.Ymean = mean(Y);
  reg.beta  = inv(X'*X)*(X'*Y);
  reg.bse   = nan(K,1);
  reg.t     = nan(K,1);
  reg.p     = nan(K,1);
  reg.bint  = nan(K,2);
  reg.df    = N-K;
  reg.alpha = 0.05;
  reg.P     = X*inv(X'*X)*X';
  reg.M     = eye(size(reg.P)) - reg.P;
  reg.e     = reg.M*Y;
  reg.Yhat  = reg.P*Y;
  reg.SSR   = reg.e'*reg.e;
  reg.s2    = reg.SSR/(N-K);
  reg.rmse  = sqrt(reg.s2);
  reg.V     = reg.s2*inv(X'*X);
  reg.Vrob  = (inv(X'*X)*X')*(reg.e*reg.e')*(X*inv(X'*X))
  reg.R2    = 1 - reg.SSR/((reg.Y-reg.Ymean)'*(reg.Y-reg.Ymean));
  reg.R2uc  = 1 - reg.SSR/(reg.Y'*reg.Y); % = (reg.Yhat'*reg.Yhat)/(reg.Y'*reg.Y)
  reg.R2adj = 1 - (1-reg.R2)*(N-1)/(N-K);

  reg.bse(:,1)  = sqrt(diag(reg.V));
  reg.t         = reg.beta./reg.bse;
  reg.p         = (1-tcdf(abs(reg.t), reg.df))*2;
  reg.bint(:,1) = reg.beta + tinv(reg.alpha/2,  reg.df)*reg.bse;
  reg.bint(:,2) = reg.beta + tinv(1-reg.alpha/2,reg.df)*reg.bse;
end

