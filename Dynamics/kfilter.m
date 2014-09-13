function [L, s_end, ss_end, s_filt, ss_filt, y_prederr, varargout] = ...
  kfilter(data, m_fwd, C, T, R, D, M, Q, s, ss, tv)
%
% Kalman Filter programming that accommodates:
% 1. Missing Data: Must be marked by a NaN in the data matrix.  if so,
%    then the corresponding row in the measurement equation will be
%    removed. 
% 2. Time-Varying Matrices: For all arrays that define the transition
%    and measurement equations, we first check if there is a third
%    dimension. If yes, then that third dimension is assumed to index
%    time. All checks are done for the matrices individually, so you
%    could have a time-varying T, but a non-time-varying M, for example.
%
%    NOTE: A matrix in the t-th location (in the third dimension of the %    array) is assumed to be the matrix that is relevant for the
%    transition equation that gets the state from time t-1 to t, or for
%    the measurement equation that relates y_t to s_t.
% 
% State Space Model Represention %
%   s_t = C + T*s_{t-1} + R*e_t (state or transition equation)
%   y_t = D + M*s_t + Q*eta_t   (observation or measurement equation)
%
% Function arguments:
%   data    (Ny x capT) matrix containing data (y(1),...,y(capT))'
%   m_fwd   The number of steps to forecast after the end of the data,
%             i.e. you forecast to capT+m
%   C       (Ns x 1) column vector representing the constant terms in
%             the transition equation
%   T       (Ns x Ns) transition matrix
%   R       (Ns x Ns) covariance matrix for exogenous shocks e_t in the
%             transition equation
%   D       (Ny x 1) vector for constant terms in the measurement
%             equation
%   M       (Ny x Ns) matrix for the measurement equation.
%   Q       (Ny x Ny) covariance matrix for the exogenous shocks eta_t
%             in the measurment equation
%   s       (Nz x 1) initial state vector.
%   ss      (Nz x Nz) covariance matrix for initial state vector
%   tv      A indicator for whether or not 1 or more matrices in the
%           state transition or measurement equation are time varying.
%           
%
% Required Output: Required because these are required for the Kalman
% Smoother, so you'll probably want them.
%
%   L           Likelihood, provided that the errors are normally
%                 distributed
%   s_end       Final filtered state vector 
%   ss_end      Covariance matrix for the final filtered state vector
%   s_filt      (Ns x t) matrix consisting of the filtered states
%   ss_filt     (Ns x Ns x t) array consisting of the filtered
%                 covariance matrices
%   y_prederr   (Ny x t) matrix consisting of the prediction error
%
%
% Optional Output: Optional because they can be used to illustrate the
% updating procedure over time. But they aren't as important as the
% required output.
%
%   y_pred      (Ny x t) matrix consisting of the y_{t+1|t}
%   yy_pred     (Ny x Ny x t) array consisting of the covariance matrix
%                 for y_{t+1|t}
%   s_pred      (Ns x t) matrix consisting of the s_{t+1|t}
%   ss_pred     (Ns x Ns x t) matrix consisting of the covariance matrix
%                 for s_{t+1|t}
%
%   This is a M-file for MATLAB.
%   - Adapted from kalcvf.m code of Iskander Karibzhanov 5-28-02.
%   - Additions include a different characterization of the state
%     transition and measurement equations, plus a way to handle
%     time-varying matrices.
%=======================================================================

  %% Basic parameters; determine sizing and loops; used often
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

  %% Check if mats/arrays are time-varying; store all mats in structure
  if tv
    mats = struct('C', C, ...
                  'T', T, ...
                  'R', R, ...
                  'D', D, ...
                  'M', M, ...
                  'Q', Q)
    matsnames = fieldnames(mats);

    % number of elements along 3rd or t dimension
    tdim_els  = arrayfun(@(mt) size(mats.(matsnames{mt}), 3), 1:length(matsnames));

    % Which indices are time varying
    tv_inds = find(tdim_els > 1);

    % Check that 3D time-varying arrays have capT elements along 3rd dim
    too_small = intersect(find(tdim_els ~= capT), tv_inds);
    if ~isempty(too_small)
      err1 = ...
        sprintf(['The following matrices have a 3rd dimension, ' ...
          'indicating time varying.\n However, they are not the ' ...
          'right size, as size(mat, 3) ~= T.\n']);
      err2 = sprintf('%s ', matsnames{too_small});
      error([err1 err2])
    end
  end

  % Pre-Allocate matrices
  s_filt    = nan(Ns, capT);
  ss_filt   = nan(Ns, Ns, capT);
  y_prederr = nan(Ny, capT);

  if nout > 6 
    y_pred  = nan(Ny, capT);
    yy_pred = nan(Ny, Ny, capT);
    s_pred  = nan(Ns, capT);
    ss_pred = nan(Ns, Ns, capT);
  end
 
  % Log likelihood
  L = 0;

  % Set up anonymous function that can do assignment of the time t
  % values. This will be useful when we want to assign T_t (the time t
  % transition matrix) to the simpler name T. We do this because T*s is
  % way more readable than T(:,:,t)*s everywhere for all matrices in the
  % procedure; plus, now we don't have to set up 3D matrices for
  % non-time-varying stuff
  asgn = @(variable,val,t) assignin('caller', variable, val(:,:,t))
  
  for t = 1:capT

    % If time varying, assign out the time time t matrix
    for mt = tv_inds
      asgn(matsnames{mt}, idx_t(mats.matsnames{mt}), t)
    end

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
    ss_pred(:,:,t) = ss;
    y_pred(:,t)    = y;
    yy_pred(:,:,t) = yy;

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
    ss_filt(:,:,t) = ss;
    
  end
  s_end = s;
  ss_end = ss;

  if nout > 6
    varargout(1) = {y_pred};
    varargout(2) = {yy_pred};
    varargout(3) = {s_pred};
    varargout(4) = {ss_pred};
  end



end  


