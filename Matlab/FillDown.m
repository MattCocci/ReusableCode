function [X] = FillDown(X, maxfill)
% FillDown - Fill values down columns, overwriting NaNs until next non-NaN
%
% X = FillDown(X) fills down numeric values, overwriting NaNs with
% the first non-NaN value directly above in the same column.
%
% X = FillDown(X, MAXFILL) will do the same, except filling down at most
% MAXFILL rows.

  ncol = size(X,2);

  % If the maxfill length not set, set it
  nrow = size(X,1);
  if ~exist('maxfill', 'var')
    maxfill = nrow;
  end

  for c = 1:ncol
    % Find non-nan indices within column c
    notnan = find(~isnan(X(:,c)));
    notnan = [notnan; nrow+1];

    % Loop over non-nan indices and fill forward
    for r = 1:(length(notnan)-1);
      start = notnan(r);
      stop  = min([notnan(r+1)-1 start+maxfill]);
      X(start:stop,c) = X(start,c);
    end
  end

end
