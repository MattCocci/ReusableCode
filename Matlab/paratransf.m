function [unconstr, constr] = paratransf(l,u)
% Function to return functions that bound/unbound parameters so you can
% run unconstrained optimization on a problem that really does have
% parameter constraints
%
% Inputs
% ------
% l  Vector, lower bounds for param vector. Inf means no lower bound
% u  Vector, upper bounds for param vector. Inf means no upper bound
%
% Outputs
% -------
% unconstr  Function that takes parameters on some constrained domain,
%           & transforms them to an unconstrained/full-real-line domain
% constr    Inverse of unconstr. Takes transformed params on
%           unconstrained domain, puts them back on constrained domain
%
% Hence theta = constr(unconstr(theta))
%

  %% Parse Bounds:

    % Given upper and lower bound vectors (u,l) for parameters, return
    % bool arrays that say which elements bounded above, below, or both
    % Entries Inf or -Inf in (l,u) denote NOT bounded below or above resp
    bd_both  = ~isinf(u) & ~isinf(l);
    bd_below = ~bd_both  & ~isinf(l);
    bd_above = ~bd_both  & ~isinf(u);


  %% Fcns to switch back and forth btw constrained/unconstrained domains

    % Convert unconstrained y in (-Inf,Inf) to x in constrained [l,u]
    to_constr_both  = @(y,l,u) 0.5*(l+u) + 0.5*((u-l).*y) ./ sqrt(1+y.^2);
    to_constr_below = @(y,l)   l + exp(y);
    to_constr_above = @(y,u)   u - exp(y);

    % Convert constrained x in [l,u] to y in unconstrained (-Inf,Inf)
    to_unconstr_both  = @(x,l,u) (2*(x-(l+u)/2) ./ (u-l)) ./ sqrt(1-(2*(x-(l+u)/2) ./ (u-l)).^2);
    to_unconstr_below = @(x,l)   log(x-l);
    to_unconstr_above = @(x,u)   log(u-x);


  %% Fcn to transform params on constrained domain to unconstrained
  function [para_unconstr] = unconstr_fcn(para)
    para_unconstr           = para;
    para_unconstr(bd_both)  = to_unconstr_both( para(bd_both),  l(bd_both), u(bd_both));
    para_unconstr(bd_below) = to_unconstr_below(para(bd_below), l(bd_below));
    para_unconstr(bd_above) = to_unconstr_above(para(bd_above), u(bd_above));
  end


  %% Inverse of above. Maps transformed params on unconstr domain to constrained
  function [para_constr] = constr_fcn(tpara)
    para_constr           = tpara;
    para_constr(bd_both)  = to_constr_both( tpara(bd_both),  l(bd_both), u(bd_both));
    para_constr(bd_below) = to_constr_below(tpara(bd_below), l(bd_below));
    para_constr(bd_above) = to_constr_above(tpara(bd_above), u(bd_above));
  end

  unconstr = @(x) unconstr_fcn(x);
  constr   = @(x) constr_fcn(x);

end
