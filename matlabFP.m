map = @(fcns, vals) cellfun(@(f) f(vals{:}), fcns);
mapc = @(fcns, vals) cellfun(@(f) f(vals{:}), fcns, 'UniformOutput', false);
iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}(); 
paren = @(x, varargin) x(varargin{:});
curly = @(x, varargin) x{varargin{:}};
