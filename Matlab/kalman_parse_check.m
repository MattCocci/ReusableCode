function [capT, Ns, Ny, C, T, R, D, M, Q] = kalman_parse_check(sysmats, data, s0, ss0)

  %% Make sure you have all of the matrices you need
  matsnames = ['C'; 'T'; 'R'; 'D'; 'M'; 'Q'];
  matsmiss  = setdiff(matsnames, cell2mat(fieldnames(sysmats)));
  if ~isempty(matsmiss)
    matsmiss  = arrayfun(@(s) [s, ','], matsmiss, 'un', 0);
    error(['Missing matrices: ', sprintf('%s', matsmiss{:})]);
  end

  % Unpack system matrices
  C = sysmats.C;
  T = sysmats.T;
  R = sysmats.R;
  D = sysmats.D;
  M = sysmats.M;
  Q = sysmats.Q;

  % Dimensions
  Ns   = size(C,1); % num states
  Ny   = size(D,1); % num observables
  capT = size(data,2); % Number of observations

  %% Check input matrix dimensions
  if size(C,2) ~= 1,     error('C must be column vector'); end
  if size(D,2) ~= 1,     error('D must be column vector'); end
  if size(data,1) ~= Ny, error('Data and D must imply the same number of observables'); end
  if size(T,1) ~= Ns,    error('T and C must imply the same number of states'); end

  if any(size(s0) ~= [Ns 1]),     error('s0 must be column vector of length Ns'); end
  if any(size(ss0,1) ~= [Ns Ns]), error('ss0 must be Ns x Ns matrix'); end

  if any([size(T,1) size(T,2)] ~= [Ns Ns]),  error('Transition matrix, T, must be square'); end
  if any([size(M,1) size(M,2)] ~= [Ny Ns]),  error('M must be Ny by Ns matrix'); end

end
