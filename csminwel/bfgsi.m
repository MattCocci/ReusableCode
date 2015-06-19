% H = bfgsi(H0,dg,dx)
%
% Description
% -----------
% This function implements the Broyden-Fletcher-Goldfarb-Shanno (BFGS)
% update of approximate inverse Hessian matrices.
%
% This is done in a quasi-Newton optimizaton scheme, and allows us to
% avoid computing the Hessian at each iteration (as we would need to in
% a true Newton scheme). Instead, we below compute approximate inverse
% Hessian matrices recursively using only gradients and the last
% approximation.
%
% Input Arguments
% ---------------
% - H0  Previous Inverse Hessian
% - dg  Previous change in gradient
% - dx  Previous change in x
%
% Updates
% -------
% 6/19/15 Add description, remove saving of matrix
% 6/8/93  Version that updates inverse hessian instead of hessian itself.
%
% Copyright by Christopher Sims 1996.  This material may be freely
% reproduced and modified.
function H = bfgsi(H0,dg,dx)

if size(dg,2) > 1, dg = dg'; end
if size(dx,2) > 1, dx = dx'; end

Hdg = H0*dg;
dgdx = dg'*dx;
if (abs(dgdx) >1e-12)
  H = H0 + (1+(dg'*Hdg)/dgdx)*(dx*dx')/dgdx - (dx*Hdg'+Hdg*dx')/dgdx;
else
  disp('bfgs update failed.')
  disp(['|dg| = ' num2str(sqrt(dg'*dg)) '|dx| = ' num2str(sqrt(dx'*dx))]);
  disp(['dg''*dx = ' num2str(dgdx)])
  disp(['|H*dg| = ' num2str(Hdg'*Hdg)])
  H=H0;
end
