function [gbar, Omegahat, J] = GMM_BuildFunctions(g_nmat, Nobs)
%
% Construct GMM building block functions from the moment conditions
%
% Inputs
% ------
% g_nmat      Function of theta s.t. g_n(theta) returns (Nobs x Nmom)
%             matrix of sample moments, i.e. returns matrix of stacked
%             g_n' for all obs n = 1,...,N
% Nobs        Number of observations in g_nmat, i.e. number of rows
%
% Output
% ------
% gbar      Function giving the mean of g_nmat across rows/obs for any
%           particular theta
% Omegahat  Function that returns the natural consistent estimator (as a
%           function of theta) for Omega=Var(g_n)=E[g_n g_n']. Formed by
%           summing the cross product (g_n g_n') across observations.
%           (Actually summing the cross prod after *demeaning* by gbar
%           at that theta)
% J         The objective function as a function of theta and weighting
%           matrix W
%

    % Function that computes gbar, the mean of the sample moment conds
    gbar = @(theta) mean(g_nmat(theta))'/Nobs;

    % Function to compute consistent est for Omega = Var(g_n(theta))
    g_nmat_demeaned = @(theta) g_nmat(theta)-repmat(gbar(theta)', Nobs, 1);
    Omegahat        = @(theta) g_nmat_demeaned(theta)'*g_nmat_demeaned(theta)/Nobs;

    % Function to compute J statistic
    J = @(theta,W) Nobs*gbar(theta)'*W*gbar(theta);
end
