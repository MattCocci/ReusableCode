% NaN out a series, x (column vec), up until we hit the nth elelemt of x
nan_before_n = @(x, n) [nan(n-1, 1); x(n:end)];

% Generate a moving average from x (column vec)
ma_gen = @(x, pds) nan_before_n(filter(ones(1, pds)/pds, 1, x), pds);

% Generate a matrix of lags
% d is a column vector; generate all lags up to n lags
% From there, you can subset what's returned to get particular lags  
lag_mat = @(d, nlags) [nan(1,nlags); toeplitz(d(1:end-1), [d(1), nan(1, nlags-1)])]
