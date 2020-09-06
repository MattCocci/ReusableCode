module OptimConstrained
#="""=#
#=Module to do *constrained* optimization using the *unconstrained*=#
#=optimization package Optim.jl. Key function exported is=#
#=optimize_constrained(). Everything else is supporting functions.=#

#=This module exports a function optimize_constrained( ), which takes=#
#=- Function f(x) where x is a N×1 vector constrained to lie in some=#
  #=restricted model-domain that is a subset of Rᴺ=#
#=- Starting point x0=#
#=- Vector in Rᴺ of lower bounds for x. Use -Inf if no lower constraint=#
#=- Vector in Rᴺ of upper bounds for x. Use  Inf if no upper constraint=#
#=- Any set of options that you might also pass to optimize( ) from=#
  #=package Optim.jl=#

#=It will then transfom the parameters from the constrained model domain=#
#=to instead take values in *all* of Rᴺ, then use optimize( ) from the=#
#=Optim.jl package to do *unconstrained* optmization on this enlarged=#
#=and unconstrained domain, and then transform the result back into the=#
#=constrained model domain to arrive at a constrained minimum of f(x)=#

#=Notation:=#
#=- Params on the constrained model domain are denoted by x or x0=#
#=- Params on the full unconstrained domain Rᴺ are denoted by y or y0=#
#="""=#

  import Optim: optimize
  export optimize_constrained


  # Fcn that returns boolean arrays (same dimension as u & l) that mark
  # parameters/indices as bounded above, below, or both
  function parse_bounds(l::AbstractArray, u::AbstractArray)
    bd_both  = ~isinf(u) & ~isinf(l)
    bd_below = ~bd_both  & ~isinf(l)
    bd_above = ~bd_both  & ~isinf(u)

    return bd_both, bd_below, bd_above
  end

  # Fcns to convert y ∈ [-∞,∞] to x ∈ [l,u]. Args are vectors
  to_model_both(y,l,u) = 0.5*(l+u) + 0.5*(u-l).*y ./ sqrt(1+y.^2)
  to_model_below(y,l)  = l + exp(y)
  to_model_above(y,u)  = u - exp(y)

  # Fcns to convert x ∈ [l,u] to y ∈ [-∞,∞]. Args are vectors
  to_unconstr_both(x,l,u) = (2*(x-(l+u)/2) ./ (u-l)) ./ sqrt(1-(2*(x-(l+u)/2) ./ (u-l)).^2)
  to_unconstr_below(x,l)  = log(x-l)
  to_unconstr_above(x,u)  = log(u-x)


  # Function for transform parameters to unconstrained domain or to
  # model domain
  function paratransf(para::AbstractArray, l::AbstractArray, u::AbstractArray,
                      bd_both::Array{Bool}, bd_below::Array{Bool}, bd_above::Array{Bool},
                      to_unconstr::Bool)

    tpara = copy(para)

    # Transform to unconstrained domain
    if to_unconstr
      tpara[bd_both]  = to_unconstr_both( para[bd_both],  l[bd_both],  u[bd_both])
      tpara[bd_below] = to_unconstr_below(para[bd_below], l[bd_below])
      tpara[bd_above] = to_unconstr_above(para[bd_above], u[bd_above])

    # Transform to constrained model domain
    else
      tpara[bd_both]  = to_model_both( para[bd_both],  l[bd_both],  u[bd_both])
      tpara[bd_below] = to_model_below(para[bd_below], l[bd_below])
      tpara[bd_above] = to_model_above(para[bd_above], u[bd_above])
    end

    return tpara
  end


  # Optimize fcn starting at x0, staying within bounds [u,l]. Can take
  # some of the optional arguments that optimize( ) does
  function optimize_constrained(fcn::Function, x0::AbstractArray,
                                l::AbstractArray, u::AbstractArray;
                                iterations=1000)

    # Figure out which parameters are bounded above, below, or both
    bd_both, bd_below, bd_above = parse_bounds(l,u)

    # Define functions for converting params back and forth
    model_to_unconstr(x) = paratransf(x, l, u, bd_both, bd_below, bd_above, true)
    unconstr_to_model(y) = paratransf(y, l, u, bd_both, bd_below, bd_above, false)

    # Translate starting point in model domain into a point of the
    # unconstrained domain
    y0 = model_to_unconstr(x0)

    # Define a closure: This function takes y on unconstrained domain,
    # transforms back to some point in the domain of the model y -> x,
    # and evaluates fcn at that point x. We will do unconstrained
    # optimization on objfcn
    objfcn(y) = fcn(unconstr_to_model(y))

    # Do unconstrained optimization on objfcn
    res = optimize(objfcn, y0, iterations=iterations)

    # Convert the values back to the model domain
    res.initial_x = unconstr_to_model(res.initial_x)
    res.minimizer = unconstr_to_model(res.minimizer)

    return res

  end

end
