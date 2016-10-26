function [est] = GMM_EstimateLinearModel(y,X,Z,twostage,W0)
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
% analytical formulas for the standard errors
%
%   betahat(W) = (X'*Z*W*Z'*X) \ (X'*Z*W*Z'y)
%
% Inputs
% ------
% y         (Nobs x 1) vector of dependent variable observations
% X         (Nobs x NXvars) matrix of (possibly endogenous) regressors
% Z         (Nobs x NZvars) matrix of predetermined instruments. Note
%           that we must have NZvars >= NXvars for identification
% twostage  Optional arg. Do two-stage efficient GMM?  Default: true
% W0        Optional arg. Initial weighting matrix. Default is eye(NZvars)
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

    % Given any matrix A, compute it's "square/cross product" A'*A
    crossprod = @(A) A'*A;

    % Function to return GMM est of params given weighting matrix W
    % NOTE: beta_init is an unused argument, but necessary for passing
    % this as a function to GMM_ImplementTwoStage
    estparams = @(W, beta_init) (X'*Z*W*Z'*X)\(X'*Z*W*Z'*y);

    % Function to compute errors given a guess for beta
    errors = @(beta) (y-X*beta);

    % Given parameter estimates beta, return constistent estimate for
    % Omega = Var(g_n(beta)) at true beta
    Omegahat = @(beta) crossprod( repmat(errors(beta),1,NZvars).*Z )/Nobs;

    % Consistent estimator for G=E[dg_n/dbeta'] always the same,
    % regardless of beta
    Ghat = @(beta) Z'*X/Nobs;

    % Function to compute J statistic
    gbar = @(beta)   Z'*(y-X*beta)/Nobs;
    J    = @(beta,W) Nobs*gbar(beta)'*W*gbar(beta);

  %% Estimate the model

    est = GMM_ImplementTwoStage(W0, NaN, estparams, J, Ghat, Omegahat, twostage);

end

