% [g, badg] = numgrad(fcn,x,varargin)
%
% Description
% -----------
% This function calculates the numerical gradient of a function.
%
% Input Arguments
% ---------------
% - fcn       A function whose numerical gradient we want to compute;
%             first argument must be the vector or scalar of parameters
%             that we would like to compute partials derivatives
%             relative to.
% - x         The point at which we want the gradient.
% - varargin  Additional arguments for fcn after x
%
% Updates
% -------
% 6/19/15 Remove commented out code; rename some things

function [g, badg] = numgrad(fcn,x,varargin)

delta = 1e-6;
Nx    = length(x);
tvec  = delta*eye(Nx);
g     = zeros(Nx,1);

f0 = feval(fcn, x, varargin{:});

badg=0;
for n = 1:Nx
  if size(x,1)>size(x,2)
    tvecv=tvec(:,n);
  else
    tvecv=tvec(n,:);
  end

  g0 = (feval(fcn, x+tvecv, varargin{:}) - f0) / delta;
  if abs(g0)< 1e15
    g(n)=g0;
  else
    disp('bad gradient ------------------------')
    g(n)=0;
    badg=1;
  end
end
