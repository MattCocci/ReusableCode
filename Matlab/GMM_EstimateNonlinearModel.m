function [est] = GMM_EstimateNonlinearModel(g_nmat, Ghat, theta_init, twostage, W0, l, u)
%
% Estimate nonlinear model via GMM. Assume model defined by moment conds
%
%   0 = E[g_n(theta)]
%
% g_n is (Nmom x 1) vector of moments that is a fcn of (Npara x 1)
% vector of parameters, denoted theta.
%
% Inputs
% ------
% g_nmat      Function of theta s.t. g_n(theta) returns (Nobs x Nmom)
%             matrix of sample moments, i.e. returns matrix of stacked
%             g_n' for all obs n=1,...,N
% Ghat        Function of theta s.t. Ghat(theta) returns sample mean of
%             dg_n/dtheta'. Necessary to construct standard errors
% theta_init  Initial guess for param vector. Starting pt for optimization
% twostage    Optional arg. Do two-stage efficient GMM? Default: true
% W0          Optional arg. Init weighting matrix. Default: eye(Nmom)
% l           Optional arg. (Npara x 1) vector of lower bounds for parameter
%             vector theta. Allows estimation with constraints.
%             Set l(i)=-Inf if ith parameter not bounded below
% u           Optional arg. (Npara x 1) Vector of upper bounds for
%             parameter vector theta. Set u(i)=Inf if ith param not
%             bounded above.


  %% Dimensions and Defaults for optional arguments

    [Nobs, Nmom] = size(g_nmat(theta_init));
    Npara        = length(theta_init);

    if ~exist('twostage', 'var')
      twostage = true;
    end
    if ~exist('W0', 'var')
      W0 = eye(Nmom);
    end
    if ~exist('l', 'var')
      l = -Inf*ones(Npara,1);
      u =  Inf*ones(Npara,1);
    end

    % Options/Defaults for fminsearch
    %options = optimset('Display','iter','PlotFcns',@optimplotfval);
    %options = optimset('Display','iter', 'MaxIter', Inf, 'MaxFunEvals', Inf, 'TolFun', 10e-15, 'TolX', 10e-15);
    options = optimset('MaxIter', Inf, 'MaxFunEvals', Inf, 'TolFun', 1e-8, 'TolX', 1e-8);


  %% Construct functions for GMM estimation

    % Function that computes gbar, the mean of the sample moment conds
    gbar = @(theta) mean(g_nmat(theta))'/Nobs;

    % Function to compute consistent est for Omega = Var(g_n(theta))
    g_nmat_demeaned = @(theta) g_nmat(theta)-repmat(gbar(theta)', Nobs, 1);
    Omegahat        = @(theta) g_nmat_demeaned(theta)'*g_nmat_demeaned(theta)/Nobs;

    % Function to compute J statistic
    J = @(theta,W) Nobs*gbar(theta)'*W*gbar(theta);


  %% Model Estimation

    % Get functions to transform parameters btw bounded/unbounded domain
    [unconstr, constr] = paratransf(l,u);

    % Estimation function
    % - Take weighting matrix + init param guess on constrained domain
    % - Convert init param guess to unconstrained domain with unconstr()
    % - fminsearch to search over unconstrained domain to min J stat
    % - Given result from fminsearch, use constr() to put the result
    %   back on constr domain
    estparams = @(W,init) constr( fminsearch(@(x) J(constr(x),W), unconstr(init), options) );


    % Estimate
    est = GMM_ImplementTwoStage(W0, theta_init, estparams, J, Ghat, Omegahat, twostage);


end


