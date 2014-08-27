function [L,zend,Pend,varargout] = ...
  kfilter(data, lead, a, F, b, H, R, Q, M, G, W, s, P)
% This version of the Kalman Filter is designed for a model with time-varying
% volatilities of the shocks. It is assumed that this filtering is being
% done conditional on a set of time-varying volatility states that "adjust" % the variance of the shocks in each period. This is done by adjusting the
% shock covariance matrix in each period.

% For this reason, additional inputs are needed. First the volatility
% states must be included, in the "W" matrix, which should be Nshocks x T. 
% Moreover, I am separating out the components of what would normally be 
% called "var" or "VVall" to make it easier to adjust the shock covariance matrix.
% Therefore R and Q will be included instead of the standard "V" matrix
% that accounts for the effects of both. Note also that this program will
% only work if the shocks are uncorrelated, i.e. Q is a diagonal matrix.

% More generally, this program is adapted from kalcvf2NaN, which is designed 
% to deal with missing data, which MUST correspond to NaN in the 'data' matrix
% if an element of the vector y(t) is missing (NaN) for the observation t, the corresponding row is ditched from the 
% measurement equation.
%
% DG 4/9/10
%
% State Space Model and Components
% 
%     z(t+1) = a+F*z(t)+R*eta(t)     (state or transition equation)
%       y(t) = b+H*z(t)+eps(t)     (observation or measurement equation)
%
%   s_t = C + T*s_{t-1} + R*e_t (state or transition equation)
%   y_t = D + M*s_t + Q*eta_t   (observation or measurement equation)
%
%   [logl, <pred, vpred, <filt, vfilt>>]=kalcvf(data, lead, a, F, b, H, var, <z0, vz0>)
%   computes the one-step prediction and the filtered estimate, as well as their covariance matrices.
%   The function uses forward recursions, and you can also use it to obtain k-step estimates.
%
%   The inputs to the KALCVF function are as follows:
%     data is a Ny?T matrix containing data (y(1), ... , y(T)).
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
    not_nan = ~isnan(data(t,:));
    Ny_t = length(data_t);
    
    data_t = data(t,not_nan);
    M_t = M(not_nan,:); 
    Q_t = Q(not_nan,not_nan);
    D_t = D(not_nan);

    %% From Filtered to Forecast values
    s  = C + T*s;           % mu_{t|t} -> mu_{t+1|t}
    ss = T*ss*T' + R;       % Sigma_{t|t} -> Sigma_{t+1|t}
    y  = D_t + M_t*s;       % E_t[y_{t+1}]
    yy = M_t*ss*M_t' + Q_t; % Var_t[y_{t+1}]

    %% Update: From Forecast to Filtered
    err = data_t - y;
    Mss = M_t*ss;
    s  = s + Mss'*(yy\err);
    ss = ss - Mss*(yy\Mss); 

    if nout > 3
      pred(:,t) = z;
      vpred(:,:,t) = P;
      if nout> 7
        yprederror(notis_nan,t) = dy;
        ystdprederror(notis_nan,t) = dy./sqrt(diag(D));
      end
    end

%     if det(D) < 10^(-4)
%         keyboard;
%     end
    ddy = D\dy;
    L = L-.5*log(det(D))-.5*dy'*ddy-.5*Ny_t*log(2*pi);
    
    
    %% updating
    PHG = (P*H_t'+G_t);
    z = z+PHG*ddy;
    P = P-PHG/D*PHG';


    if nout > 5
      filt(:,t) = z;
      vfilt(:,:,t) = P;
    end
  end
  zend = z;
  Pend = P;
  
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


