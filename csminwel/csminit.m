function [fhat,xhat,fcount,retcode] = csminit(fcn,x0,f0,g0,badg,H0,varargin)
% CSMINIT - Performs Berndt, Hall, Hall, & Hausman 1974 (BHHH) line
%           search to find lambda, the optimal step size
%
% This function will take a starting point x0, along with the value of
% the function, the gradient, and the approximate inverse Hessian H0
% (all at x0) to determine a new point to move to.
%
% This function determines the point to move to by using the approximate
% inverse Hessian times the gradient to get a Quasi-Newton step
% direction. If this direction is insufficiently aligned with the
% gradient, a "low angle correction" is introduced, as in BHHH. (This is
% the "Restriction on Q" section.)
%
% After determining the step direction, perform a line search to
% determine a good value of lambda such that f(x0 + lambda*direction)
% yields sufficient improvement. (This is the "Criterion for choice of
% lambda" section in BHHH)
%
% retcodes document the results of the line search.
%
% INPUT ARGUMENTS
% ---------------
% fcn       Function to be minimized
% x0        Point we're moving from
% f0        Value of function at point we're moving from
% g0        Gradient at point we're moving from
% badg      Indicator for bad gradient
% H0        Approximation for the inverse hessian at the point we're
%           moving from. In combination with the gradient, this will be
%           used to determine the step direction
% varargin  Additional inputs to fcn aside from parameters (1st arg)
%
%
% retcodes
% --------
% 0   normal step
% 1   zero gradient
% 2,4 back and forth adjustment of stepsize didn't finish
% 3   smallest stepsize still improves too slow.
% 5   largest step still improves too fast.
% 6   no improvement found.
%
% Changes
% -------
% Updated 6/21/15 by MDC to provide more comments about the code, remove
% commented out code, eliminate unnecessary variables, and rename things
% more suggestively to be in line with the BHHH paper that this line
% search is based on.
%
% Fixed 7/17/93 to use inverse-hessian instead of hessian itself in bfgs
% update.

% BHHH settings
alpha = .005; % previously named ANGLE; now corresponds to alpha in BHHH
delta = .3;   % previously named THETA; now corresponds to delta in BHHH;
              %   (0 < delta <.5) delta near .5 makes long line searches,
              %   possibly fewer iterations.

FCHANGE = 1000;
MINLAMB = 1e-9;
MINDFAC = .01;

% Size of the gradient at x0
gnorm = norm(g0);

% Pre-define variables that will be returned
fcount    = 0;  % Number of times fcn is evaluated through the line search
xhat      = x0; % Will be the new suggested x value
fhat      = f0; % Will be the fcn at the new suggested x value
lambdahat = 0;  % Will be the accepted lambda that gives xhat and fhat

if (gnorm < 1.e-12) & ~badg % gradient convergence
  retcode =1;
  dxnorm=0;
else

  % Get the Quasi-Newton recommended direction, dx, and the predicted
  % improvement (Newton direction times gradient)
  dx = -H0*g0;
  dxnorm = norm(dx);
  if dxnorm > 1e12
    disp('Near-singular H problem.')
    dx = dx*FCHANGE/dxnorm;
  end
  dfhat = dx'*g0;


  % Test for alignment of dx with gradient and fix if necessary
  if ~badg
    a = -dfhat/(gnorm*dxnorm);
    if a < alpha
      dx = dx - (alpha*dxnorm/gnorm+dfhat/(gnorm*gnorm))*g0;
      dx = dx*dxnorm/norm(dx)    % This keeps scale invariant to the angle correction
      dfhat = dx'*g0;
      disp(sprintf('Correct for low angle: %g',a))
    end
  end
  disp(sprintf('Predicted improvement: %18.9f',-dfhat/2))

  % Have OK dx, now adjust length of step (lambda) until min and max
  % improvement rate criteria are met. (Do the line search)
  done       = 0;
  factor     = 3;
  shrink     = 1;
  lambda     = 1;
  lambdaMin  = 0;
  lambdaMax  = Inf;
  lambdaPeak = 0;
  fPeak      = f0;
  while ~done
    if size(x0,2)>1
      dxtest=x0+dx'*lambda;
    else
      dxtest=x0+dx*lambda;
    end

    % Get the current value of the function
    f = feval(fcn,dxtest,varargin{:});
    disp(sprintf('lambda = %10.5g; f = %20.7f',lambda,f ))

    % If you have an improvement, update the variables to be returned to
    % catch this
    if f<fhat
      fhat      = f;
      xhat      = dxtest;
      lambdahat = lambda;
    end
    fcount=fcount+1;

    % Setting lambda:
    % Shrink if no or insufficient (< delta, see BHHH) improvement
    % Grow if improvement > 1-delta (see BHHH)
    shrinkSignal = (~badg & (f0-f < max([-delta*dfhat*lambda 0]))) | (badg & (f0-f) < 0) ;
    growSignal   = ~badg & ( (lambda > 0)  &  (f0-f > -(1-delta)*dfhat*lambda) );

    if shrinkSignal & ( (lambda > lambdaPeak) | (lambda < 0) )
      if (lambda > 0) & ((~shrink) | (lambda/factor <= lambdaPeak))
        shrink=1;
        factor=factor^.6;
        while lambda/factor <= lambdaPeak
          factor=factor^.6;
        end
        if abs(factor-1) < MINDFAC
          if abs(lambda) < 4
            retcode=2;
          else
            retcode=7;
          end
          done=1;
        end
      end
      if (lambda < lambdaMax) & (lambda > lambdaPeak)
        lambdaMax=lambda;
      end
      lambda = lambda/factor;
      if abs(lambda) < MINLAMB
        if (lambda > 0) & (f0 <= fhat)
          % try going against gradient, which may be inaccurate
          lambda = -lambda*factor^6
        else
          if lambda < 0
            retcode = 6;
          else
            retcode = 3;
          end
          done = 1;
        end
      end
    elseif (growSignal & lambda>0) | (shrinkSignal & ((lambda <= lambdaPeak) & (lambda > 0)))
      if shrink
        shrink=0;
        factor = factor^.6;
        if abs(factor-1)<MINDFAC
          if abs(lambda)<4
            retcode=4;
          else
            retcode=7;
          end
          done=1;
        end
      end
      if (f < fPeak) & (lambda>0)
        fPeak=f;
        lambdaPeak=lambda;
        if lambdaMax <= lambdaPeak
          lambdaMax=lambdaPeak*factor*factor;
        end
      end
      lambda=lambda*factor;
      if abs(lambda) > 1e20;
        retcode = 5;
        done =1;
      end
    else
      done=1;
      if factor < 1.2
        retcode=7;
      else
        retcode=0;
      end
    end
  end
end
disp(sprintf('Norm of dx %10.5g', dxnorm))

