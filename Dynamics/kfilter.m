function [L,zend,Pend,varargout] = ...
  kfilter(data, m_fwd, C, T, R, D, M, Q, s, P, tv)
%
% Kalman Filter programming that accommodates:
% 1. Missing Data: Must be marked by a NaN in the data matrix.
%    if so, then the corresponding row in the measurement
%    equation will be removed. 
% 2. Time-Varying Matrices: For all arrays that define the
%    transition and measurement equations, we first check if
%    there is a third dimension. If yes, then that third
%    dimension is assumed to index time. All checks are done
%    for the matrices individually, so you could have a
%    time-varying T, but a non-time-varying M, for example.
%
%    NOTE: A matrix in the t-th location (in the third dimension
%    of the array) is assumed to be the matrix that is relevant
%    for the transition equation that gets the state from time
%    t-1 to t, or for the measurment equation that relates y_t
%    to s_t.
% 
% State Space Model Represention
%
%   s_t = C + T*s_{t-1} + R*e_t (state or transition equation)
%   y_t = D + M*s_t + Q*eta_t   (observation or measurement equation)
%
% Function arguments:
%   data  (Ny x capT) matrix containing data (y(1),...,y(capT))'
%   m_fwd The number of steps to forecast after the end of the
%            data, i.e. you forecast to capT+m
%   C     (Ns x 1) column vector representing the constant terms
%           in the transition equation
%   T     (Ns x Ns) transition matrix
%   R     (Ns x Ns) covariance matrix for exogenous shocks e_t
%           in the transition equation
%   D     (Ny x 1) vector for constant terms in the measurement
%           equation
%   M     (Ny x Ns) matrix for the measurement equation.
%   Q     (Ny x Ny) covariance matrix for the exogenous shocks
%           eta_t in the measurment equation
%   z is an optional Nz?1 initial state vector.
%   P is an optional Nz?Nz covariance matrix of an initial state vector.
%   tv    A indicator for whether or not 1 or more matrices in
%           the state transition or measurement equation are
%           time varying. If not, save some computational costs
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
%   Adapted from code of Iskander Karibzhanov 5-28-02.
%==========================================================================

  capT = size(data,2);
  Ns   = size(C,1);
  Ny   = size(D,1);

  nout = nargout;

  %% Check input matrix dimensions
  if size(C,2) ~= 1, error('C must be column vector') end
  if size(D,2) ~= 1, error('D must be column vector') end
  if size(data,2) ~= Ny
    error('Data and D must have the same number of rows')
  end
  if any([size(T,1), size(T,2)] ~= [Ns Ns])
    error('Transition matrix, T, must be square')
  end
  if any([size(M,1) size(M,2)] ~= [Ny Ns])
    error('M must be Ny by Ns matrix')
  end
  if size(T,1) ~= Ns
    error('T and C must have the same number of rows')
  end
  if any(size(s) ~= [Ns 1])
    error('s0 must be column vector of length Ns')
  end
  if any(size(ss,1) ~= [Ns Ns])
    error('ss0 must be Ns by Ns matrix') end
  if any(R ~= diag(diag(R)))
      error('R must be a diagonal matrix')
  end
  if any(Q ~= diag(diag(Q)))
      error('Q must be a diagonal matrix')
  end

  %% Check which matrices/arrays are time-varying
  if tv
    matsnames = {'C', 'T', 'R', 'D', 'M', 'Q'};
    matscell  = {C, T, R, D, M, Q};
    tv_inds = find(cellfun(@(m) ndims(m) > 2, mats, ...
                   'UniformOuput', true));
    for mt = 1:length(matsnames)
      mats.(matsnames{mt}) = matscell{mt};
    end
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
 
  % Log likelihood
  L = 0;
  
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

    %% Evaluate the likelihood p(y_t | I_{t-1},C,T,Q,D,M,R)
    err = data_t - y;     
    L = L - 0.5*Ny_t*log(2*pi) - 0.5*log(det(yy)) ...
        - 0.5*err'*(yy\err);

    %% From Forecast to Filtered
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


