%% numgrad.m
%
% Compute numerical gradient of fcn at x.
%
% Updated version of numgrad.m by Chris Sims and Jinill Kim.
function [g, badg] = numgrad(fcn, x, args, opt)

% Default settings
defaults = struct(...
            'delta',     1e-6, ... % Size of perturbations
            'two_sided', 1, ...
            );

% For those settings not provided, use defaults
if ~exist('var', 'opt')
  opt = defaults;
else
  flds = fieldnames(defaults);
  for f = 1:length(flds)
    if ~isfield(opt, flds{f})
      opt.(flds{f}) = defaults.(flds{f});
    end
  end
end


n    = length(x);
tvec = delta*eye(n);

f0 = feval(fcn, x, varargin{:});

g = zeros(n,1);

if

for i=1:n
   if size(x,1)>size(x,2)
      tvecv=tvec(i,:);
   else
      tvecv=tvec(:,i);
   end

   g0 = (feval(fcn, x+tvecv, varargin{:}) - f0) ...
         /(delta);
   if abs(g0)< 1e15
      g(i)=g0;
   else
      disp('bad gradient ------------------------')
      g(i)=0;
      badg=1;
   end
end
