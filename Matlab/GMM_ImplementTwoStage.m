function [est] = GMM_ImplementTwoStage(N, W0, theta_init, estparams, J, Ghat, Omegahat, twostage)
%
% Implment two stage GMM procedure to obtain param estimates, standard
% errors, J statistic. Returns structure array with objects, first entry
% the first stage estimates, and second entry the second stage
% estimates.
%
% Inputs
% ------
% W0          Initial weigting matrix
%
% theta_init  Initial guess for parameter vector
%
% estparams   Function s.t. estparams(W,theta_init) returns an estimate
%             for parameter vector theta.
%
%             If called by GMM_EstimateNonlinearModel, this function is
%             a wrapper function for a numerical optimization routine
%             which uses theta_init as starting point for numerical
%             search. If called by GMM_EstimateLinearModel, the
%             estimator is available in closed form, hence theta_init is
%             a nuisance argument and estparams(W,theta_init) directly
%             computes the estimate by formula rather than numerical
%             optimization
%
% J           Function s.t. J(theta,W) returns the value of the
%             objective function aka the J-statistic at theta for
%             weighted matrix W
%
% Ghat        Function s.t. Ghat(theta) returns consistent sample
%             estimate for the expected Jacobian matrix of moment
%             conditions, i.e. sample analog of E[dg_n(theta)/dtheta'].
%             Needed for constructing standard errors
%
% Omegahat    Function s.t. Omegahat(theta) gives consistent sample
%             estimate of the covariance matrix of moment conditions,
%             i.e. sample analog of Var(g_n(theta))
%
% twostage    Whether to also report the two-step estimate using
%             Omegahat^-1 at the first step estimate as the weighting
%             matrix
%

    % Initialize strutcure
    struct(1+twostage).thetahat = [];

    % Compute (and store) first stage parameter estimates and stderrs
    est(1).thetahat = estparams(W0, theta_init);
    G               = Ghat(est(1).thetahat);
    Omega           = Omegahat(est(1).thetahat);
    est(1).avar     = (G'*W0*G)\(G'*W0*Omega*W0*G)*inv(G'*W0*G);
    est(1).stderrs  = sqrt(diag(est(1).avar)/N);
    est(1).J        = J(est(1).thetahat, W0);
    est(1).Ghat     = G;
    est(1).Omegahat = Omega;

    % If also want estimates and standard errors for two-stage estimate
    if twostage
      W1              = inv(Omega); % Optimal weighting matrix
      est(2).thetahat = estparams(W1, est.thetahat0);
      G               = Ghat(est(2).thetahat);
      Omega           = Omegahat(est(2).thetahat);
      est(2).avar     = inv(G'*(Omega\G));
      est(2).stderrs  = sqrt(diag(est(2).avar)/N);
      est(2).J        = J(est(2).thetahat,W1);
      est(2).Ghat     = G;
      est(2).Omegahat = Omega;
    end
end
