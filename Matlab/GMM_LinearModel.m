function [est] = GMM_LinearModel(y,X,Z,twostage,W0)
%
% Function to estimate single equation linear regression model with
% potentially endogenous regressors via GMM. Model assumed to be of form
%
%   y_n = x_n'*beta + e_n  where E[e_n*z_n] = 0
%
% where z_n is a vector of predetermined instruments orthogonal to the
% error term e_n. The identification assumption that E[e_n*z_n]=0 serves
% as a set of moment conditions that we can exploit to estimate the
% model via GMM.
%
% Specifically, define GMM estimate betahat(W) cond'l on weighting mat W
%
%   betahat(W) := min_b J(b,W)
%    J(b,W)    := Nobs*gbar(b)'*W*gbar(b)
%   gbar(b)    := Z'(y-Xb)/Nobs
%
% The FOCs of the above min problem lead to the following analytical
% formulas for the estimates, while asymptotic theory provides the
% following analytical formulas for the standard errors
%
%   betahat(W) = (X'*Z*W*Z'*X) \ (X'*Z*W*Z'y)
%
% Inputs
% ------
% y         (Nobs x 1) vector of dependent variable observations
% X         (Nobs x NXvars) matrix of potentially endogenous regressors
% Z         (Nobs x NZvars) matrix of predetermined instruments. Note
%           that we must have NZvars >= NXvars for identification
% twostage  Optional. Whether do two-stage efficient GMM. Default true
% W0        Optional. Initial weighting matrix. Default is eye(NZvars)
%
%

  if ~exist('twostage', 'var')
    twostage = true;
  end

  % Get dimensions of inputs
  [Nobs, NXvars] = size(X);
  NZvars         = size(Z,2);
  if ~exist('W0', 'var')
    W0 = eye(NZvars);
  end

  %% General purpose functions

    % Given any matrix A, compute "square" A'*A
    squaremat = @(A) A'*A;

    % GMM estimate of parameters weighting matrix W
    estparams = @(W) (X'*Z*W*Z'*X)\(X'*Z*W*Z'*y);

    % Compute errors
    errors = @(beta) (y-X*beta);

    % Given parameter estimates beta, return constistent estimate for
    % Omega, the asymptotic var of g_n(beta)
    Omegahat = @(beta) squaremat( repmat(errors(beta),1,NZvars).*Z )/Nobs;

    % Consistent estimator for G=E[dg_n/dbeta'] always the same
    Ghat = Z'*X/Nobs;

    % Compute the J statistic
    gbar = @(beta)   Z'*(y-X*beta)/Nobs;
    J    = @(beta,W) Nobs*gbar(beta)'*W*gbar(beta);

  %% Estimate the model

    % Compute (and store) first stage parameter estimates and stderrs
    est.betahat0  = estparams(W0);
    Omegahat0     = Omegahat(est.betahat0);
    est.stderrs0  = (Ghat'*W0*Ghat)\(Ghat'*W0*Omegahat0*W0*Ghat)*inv(Ghat'*W0*Ghat);
    est.J0        = J(est.betahat0,W0);

    % If also want estimates and standard errors for two-stage estimate
    if twostage
      W1            = inv(Omegahat0);
      est.betahat1  = estparams(W1);
      Omegahat1     = Omegahat(est.betahat1);
      est.stderrs1  = inv(Ghat'*W1*Ghat);
      est.J1        = J(est.betahat1,W1);
    end

end

