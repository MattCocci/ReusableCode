%% numgrad.m
%
% Compute numerical gradient of scalar-valued fcn at x, a vector of
% length n.
%
% Updated version of numgrad.m by Chris Sims and Jinill Kim.
function [g, badg] = numgrad(fcn, x, args, varargin)

% Setup
n    = length(x);
if size(x,1) < size(x,2)
  x = x';
end

% Default settings
defaults = struct(...
            'deltas',    1e-6*ones(n,1), ... % Size of perturbations
            'idx_grad',  1:n, ...
            'two_sided', 1 ...
            );

% Set up options, filling in defaults where none provided
opt = FillDefaultSettings(defaults, varargin{:});

% Amounts to perturb x by
delta = diag(opt.deltas);

% Define gradient calculating function depending upon option used
if defaults.two_sided
  grad = @(i_) (1/delta(i_,i_))*(feval(fcn, x+delta(:,i_)/2, args{:}) ...
                                  - feval(fcn, x-delta(:,i_)/2, args{:}));
else
  f0 = feval(fcn, x, args{:});
  grad = @(i_) (1/delta(i_,i_))*(feval(fcn, x+delmat(:,i_)/2, args{:}) - f0);
end

g = nan(n,1);
for i_ = opt.idx_grad

  g0 = grad(i_);

  if abs(g0)< 1e15
    g(i_)=g0;
  else
    fprintf('-- Bad gradient for parameter %d --\n', i_);
    g(i_)=0;
    badg=1;
  end
end
