function [L,zend,Pend,varargout] = kfilter(data, lead, C, T, R, D, M, Q, M, G, W, s, P)
%
% More generally, this program is adapted from kalcvf2NaN, which is designed to
% deal with missing data, which MUST correspond to NaN in the 'data' matrix if
% an element of the vector y(t) is missing (NaN) for the observation t, the
% corresponding row is ditched from the measurement equation.
%
% DG 4/9/10
%
% State Space Model and Components
%
%   s_t = C + T*s_{t-1} + R*e_t (state or transition equation)
%   y_t = D + M*s_t + Q*eta_t   (observation or measurement equation)
%
%   The inputs to the function are as follows:
%     data is a Ny x T matrix containing data (y(1), ... , y(T)).
%     lead is the number of steps to forecast after the end of the data.
%        a is an Nz?1 vector for a time-invariant input vector in the transition equation.
%        F is an Nz?Nz matrix for a time-invariant transition matrix in the transition equation.
%        b is an Ny?1 vector for a time-invariant input vector in the measurement equation.
%        H is an Ny?Nz matrix for a time-invariant measurement matrix in the measurement equation.
%        R is an Nz x Nshocks matrix that translates the exogenous shocks
%           into states in the transition equation
%        Q is the Nshocks x Nshocks covariance matrix of the exogenous
%           shocks eta_t
%        M is the Ny x Ny covariance matrix of the measurement error
%           epsilon_t
%        G is the Nz x Ny covariance matrix of R*eta_t and epsilon_t
%        W is the Nz x T matrix of time-varying volatility adjustments
%        z is an optional Nz?1 initial state vector.
%        P is an optional Nz?Nz covariance matrix of an initial state vector.
%
%   The KALCVF function returns the following output:
%     logl is a value of the average log likelihood function of the SSM
%             under assumption that observation noise eps(t) is normally distributed
%     pred is an optional Nz?(T+lead) matrix containing one-step predicted state vectors.
%    vpred is an optional Nz?Nz?(T+lead) matrix containing mean square errors of predicted state vectors.
%     filt is an optional Nz?T matrix containing filtered state vectors.
%    vfilt is an optional Nz?Nz?T matrix containing mean square errors of filtered state vectors.
%
%   The initial state vector and its covariance matrix of the time invariant Kalman filters
%   are computed under the stationarity condition:
%          z0 = (I-F)\a
%         vz0 = (I-kron(F,F))\V(:)
%   where F and V are the time invariant transition matrix and the covariance matrix of transition equation noise,
%   and vec(V) is an Nz^2?1 column vector that is constructed by the stacking Nz columns of matrix V.
%   Note that all eigenvalues of the matrix F are inside the unit circle when the SSM is stationary.
%   When the preceding formula cannot be applied, the initial state vector estimate is set to a
%   and its covariance matrix is given by 1E6I. Optionally, you can specify initial values.
%
%   This is a M-file for MATLAB.
%   Copyright 2002-2003 Federal Reserve Bank of Atlanta
%   $Revision: 1.2 $  $Date: 2003/03/19 19:16:17 $
%   Iskander Karibzhanov 5-28-02.
%   Master of Science in Computational Finance
%   Georgia Institute of Technology
%==========================================================================
% Revision history:
%
%  03/19/2003  -  algorithm and interface were adapted from SAS/IML KALCVF subroutine for use in MATLAB M file
%
%==========================================================================

  capT = size(data,2);
  Ns   = size(C,1);
  Ny   = size(D,1);

  if nin~=13
    error('Thirteen input arguments required.')
  end

  nout = nargout;

  % Check input matrix dimensions
  if size(C,2) ~= 1, error('C must be column vector') end
  if size(D,2) ~= 1, error('D must be column vector') end
  if size(data,2) ~= Ny, 
    error('Data and D must have the same number of columns')
  end
  if any(size(T) ~= [Ns Ns])
    error('Transition matrix, T, must be square')
  end
  if any(size(M) ~= [Ny Ns])
    error('M must be Ny by Ns matrix')
  end
  if any(size(s)~=[Ns 1])
    error('s0 must be column vector of length Nz')
  end
  if any(size(P)~=[Ns Ns])
    error('s0 must be Ns by Ns matrix') end
  if any(Q ~= diag(diag(Q)))
      error('Q must be a diagonal matrix')
  end

  if nout>3
    pred = zeros(Nz,T);
    vpred = zeros(Nz,Nz,T);
    if nout > 5
        filt = zeros(Nz,T);
        vfilt = zeros(Nz,Nz,T);
        if nout > 7
            yprederror = NaN*zeros(Ny,T);
            ystdprederror = NaN*zeros(Ny,T);
        end
    end
  end
 
  L = 0;
  
  qdiag = sqrt(diag(Q));

  for t = 1:capT

    %% Handle missing observations
      
      % if an element of the vector y_t is missing (NaN) for
      % the observation t, the corresponding row is ditched
      % from the measurement equation.
      not_nan = ~isnan(data(:,t));
      Ny_t = length(data_t);
      
      data_t = data(not_nan,t);
      M_t = M(not_nan,:); 
      Q_t = Q(not_nan,not_nan);
      D_t = D(not_nan);

    %% From Filtered to Forecast values
    s  = C + T*s;           % mu_{t|t} -> mu_{t+1|t}
    ss = T*ss*T' + R;       % Sigma_{t|t} -> Sigma_{t+1|t}
    y  = D_t + M_t*s;       % E_t[y_{t+1}]
    yy = M_t*ss*M_t' + Q_t; % Var_t[y_{t+1}]

    %% Save forecasts
    s_pred(:,t)    = s;
    s_vpred(:,:,t) = ss;
    y_pred(:,t)    = y;
    y_vpred(:,:,t) = yy;

    %% Update: From Forecast to Filtered
    err = data_t - y;     
    Mss = M_t*ss;           
    s  = s + Mss'*(yy\err); % mu_{t+1|t} -> mu_{t+1|t+1} 
    ss = ss - Mss*(yy\Mss); % Sigma_{t+1|t} -> Sigma_{t+1|t+1}

    %% Save filtering information
    y_prederr(:,t) = err;
    s_filt(:,t)    = s;
    s_vfilt(:,:,t) = ss;
    
  end
  s_end = s;
  ss_end = ss;
  
  if nout > 3
    varargout(1) = {pred};
    varargout(2) = {vpred};
    if nout>5
        varargout(3) = {filt};
        varargout(4) = {vfilt};
        if nout > 7   
          varargout(5) = {yprederror};
          varargout(6) = {ystdprederror};
          varargout(7) = {sqrt(mean(yprederror.^2,2))'};
          varargout(8) = {sqrt(mean(ystdprederror.^2,2))'};
        end
    end
  end


