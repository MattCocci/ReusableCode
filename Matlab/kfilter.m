function [L, KF] = kfilter(data, sysmats, s0, ss0, fullout)
% Kalman Filter program that accommodates missing data.
%
% Missing data must be marked by a NaN in the data matrix. If so, then
% the corresponding row in the measurement equation will be removed for
% that time period.
%
% State Space Model Represention
%   y_t = D + M*s_t + Q*eta_t   (observation or measurement equation)
%   s_t = C + T*s_{t-1} + R*e_t (state or transition equation)
%
% Function arguments:
%   data      (Ny x capT) matrix containing data (y(1),...,y(capT))'
%   sysmats   Struct that holds the system matrices. Must have fields
%     C       (Ns x 1) column vector representing the constant terms in
%             the transition equation
%     T       (Ns x Ns) transition matrix
%     R       (Ns x Ns) covariance matrix for exogenous shocks e_t in
%               the transition equation
%     D       (Ny x 1) vector for constant terms in the measurement
%               equation
%     M       (Ny x Ns) matrix for the measurement equation.
%     Q       (Ny x Ny) covariance matrix for the exogenous shocks eta_t
%               in the measurment equation
%   s0      (Nz x 1) initial state vector.
%   ss0     (Nz x Nz) covariance matrix for initial state vector
%   fullout Indicator for whether to return the full gamut of output
%           (notably information about the updating procedure)
%
% Automatic Output: Just the likelihood in the event that you are
% simplying mode finding.
%   L           Likelihood, provided that the errors are normally
%                 distributed
%
% Optional Output (if nargout > 1):
%   s_end       Final filtered state vector
%   ss_end      Covariance matrix for the final filtered state vector
%   s_filt      (Ns x t) matrix consisting of the filtered states
%   ss_filt     (Ns x Ns x t) array consisting of the filtered
%                 covariance matrices
%   y_prederr   (Ny x t) matrix consisting of the prediction error %
%
% Optional Output (if nargout > 1 & fullout==1): Optional because they
% can be used to illustrate the updating procedure over time. But they
% aren't as important as the optional output above.
%
%   y_pred      (Ny x t) matrix consisting of the y_{t+1|t}
%   yy_pred     (Ny x Ny x t) array consisting of the covariance matrix
%                 for y_{t+1|t}
%   s_pred      (Ns x t) matrix consisting of the s_{t+1|t}
%   ss_pred     (Ns x Ns x t) matrix consisting of the covariance matrix
%                 for s_{t+1|t}
%
%   This is a M-file for MATLAB.
%   - This Kalman Filter code is adapted from that of the FRBNY DSGE
%     model code (September 2014, Liberty Street Economics). That, in
%     turn, built upon  kalcvf.m by Iskander Karibzhanov 5-28-02.
%   - Additions include a different characterization of the state
%     transition and measurement equations
%
%=======================================================================

  nout = nargout;

  %% Return system matrices if they pass a dimension check
  [capT, Ns, Ny, C, T, R, D, M, Q] = kalman_parse_check(sysmats, data, s0, ss0);

  % Pre-Allocate matrices
  KF.s_filt    = nan(Ns, capT);
  KF.ss_filt   = nan(Ns, Ns, capT);
  KF.y_prederr = nan(Ny, capT);

  if fullout
    KF.s_pred  = nan(Ns, capT);
    KF.ss_pred = nan(Ns, Ns, capT);
    KF.y_pred  = nan(Ny, capT);
    KF.yy_pred = nan(Ny, Ny, capT);
  end

  % Set up initial state vector and cov matrix
  s  = s0;
  ss = ss0;

  % Log likelihood
  L = 0.0;

  for t = 1:capT

    %% Handle missing observations

      % if an element of the vector y_t is missing (NaN) for
      % the observation t, the corresponding row is ditched
      % from the measurement equation.
      not_nan = ~isnan(data(:,t));

      data_t  = data(not_nan,t);
      Ny_t    = length(data_t);
      M_t     = M(not_nan,:);
      Q_t     = Q(not_nan,not_nan);
      D_t     = D(not_nan);

    %% From Filtered to Forecast values
    s  = C + T*s;           % mu_{t-1|t-1} -> mu_{t|t-1}
    ss = T*ss*T' + R;       % Sigma_{t-1|t-1} -> Sigma_{t|t-1}
    y  = D_t + M_t*s;       % y_{t|t-1} = E_{t-1}[y_t]
    yy = M_t*ss*M_t' + Q_t; % Var_{t-1}[y_t]

    %% Save forecasts
    if fullout
      KF.s_pred(:,t)    = s;
      KF.ss_pred(:,:,t) = ss;
      KF.y_pred(:,t)    = D + M*s;     % <- Compute forecasts for all obs,
      KF.yy_pred(:,:,t) = M*ss*M' + Q; % <- not just the non-missing ones
    end

    %% Evaluate the likelihood p(y_t | I_{t-1},C,T,Q,D,M,R)
    err = data_t - y;
    L = L - 0.5*Ny_t*log(2*pi) - 0.5*log(det(yy)) ...
          - 0.5*err'*(yy\err);

    %% From Forecast to Filtered
    Mss = M_t*ss;
    s   = s + Mss'*(yy\err);  % mu_{t|t-1} -> mu_{t|t}
    ss  = ss - Mss'*(yy\Mss); % Sigma_{t|t-1} -> Sigma_{t|t}

    %% Save filtering information
    KF.y_prederr(:,t) = err;
    KF.s_filt(:,t)    = s;
    KF.ss_filt(:,:,t) = ss;

  end

  KF.L      = L;
  KF.s_end  = s;
  KF.ss_end = ss;


end


