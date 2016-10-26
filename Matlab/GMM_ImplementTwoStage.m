function [est] = GMM_ImplementTwoStage(W0, theta_init, estparams, J, Ghat, Omegahat, twostage)
%
% Implment two stage GMM procedure to obtain param estimates, standard
% errors, J statistic
%
% Inputs
% ------
% W0          Initial weigting matrix
% theta_init  Initial guess for parameter vector
% estparams   Function s.t. estparams(W,theta_init) returns an estimate
%             for parameter vector theta. If called by
%             GMM_EstimateNonlinearModel, this function is a wrapper
%             function for a numerical optimization routine which uses
%             theta_init as starting point for numerical search
% J           Function s.t. J(theta,W) returns the J-statistic
% Ghat        Function s.t. Ghat(theta) returns consistent sample
%             estimate for expected Jacobian matrix of moment
%             conditions, i.e. sample analog of E[dg_n(theta)/dtheta'].
%             Needed for constructing standard errors
% Omegahat    Function s.t. Omegahat(theta) gives consistent sample
%             estimate of the covariance matrix of moment conditions,
%             i.e. sample analog of Var(g_n(theta))
% twostage    Whether really to do two-stages or just stop at one stage

    % Compute (and store) first stage parameter estimates and stderrs
    est.thetahat0 = estparams(W0, theta_init);
    Ghat0         = Ghat(est.thetahat0);
    Omegahat0     = Omegahat(est.thetahat0);
    est.avar0     = (Ghat0'*W0*Ghat0)\(Ghat0'*W0*Omegahat0*W0*Ghat0)*inv(Ghat0'*W0*Ghat0);
    est.stderrs0  = sqrt(diag(est.avar0));
    est.J0        = J(est.thetahat0, W0);

    % If also want estimates and standard errors for two-stage estimate
    if twostage
      W1            = inv(Omegahat0);
      est.thetahat1 = estparams(W1, est.thetahat0);
      Ghat1         = Ghat(est.thetahat1);
      Omegahat1     = Omegahat(est.thetahat1);
      est.avar1     = inv(Ghat1'*inv(Omegahat1)*Ghat1);
      est.stderrs1  = sqrt(diag(est.avar1));
      est.J1        = J(est.thetahat1,W1);
    end
end
