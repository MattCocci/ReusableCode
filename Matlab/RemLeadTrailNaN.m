function [X] = RemLeadTrailNaN(X, choice)

  nrows = size(X, 1);
  allnan_rows = (sum(~isnan(X),2) == 0);

  if ~any(strcmp({'all', 'trail', 'lead'}, choice))
    error('Invalid option: Please select ''lead'', ''trail'', or ''all''')
  end

  if any(strcmp({'all', 'trail'}, choice))
    trailing_all_nan_rows = ((1:nrows)' == cumsum(allnan_rows(end:-1:1)));
    trailing_all_nan_rows = trailing_all_nan_rows(end:-1:1);
    X(trailing_all_nan_rows,:) = [];
  end

  if any(strcmp({'all', 'lead'}, choice))
    leading_all_nan_rows = ((1:nrows)' == cumsum(allnan_rows));
    X(leading_all_nan_rows,:) = [];
  end

end
