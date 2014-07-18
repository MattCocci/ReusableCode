
% NaN out a series, x, up until we hit the nth elelemt of x
nan_before_n = @(x, n) [nan(n-1, 1); x(n:end)];

% Generate a moving average process 
ma_gen = @(col_vec, pds) nan_before_n(filter(ones(1, pds)/pds, 1, col_vec), pds);

% Generate a matrix of lags; d is a column vector, n is how far to lag it
lag_mat = @(d, nlags) [nan(1,nlags); toeplitz(d(1:end-1), [d(1), nan(1, nlags-1)])]
