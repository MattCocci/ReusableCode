import sympy as sp
import os

# Initialize symbols
X, mu, sigma = sp.symbols('X mu sigma')

# Geometric Brownian Motion
#b = mu*X
#s = sigma*X
#f = sp.ln(X)

# Some Other Example
b = 0
s = 1
f = X**2


def invertInteractive(f, X, Y):
  """
  Given function Y = f(X) return inverse f^{-1}(X).

  Note:
  - f could be fcn of multiple variables. This inverts in variable X
  - f must be a sympy function of sympy symbol/variable X
  - Returns function of Y and anything else f depends upon, removing
    dependence upon X
  - Accounts for non-unique inverse by prompting user to select
    from possible inverses
  """

  # Try to construct inverse
  finv = sp.solve(Y-f, X)

  # If inverse unique, return it. If not, prompt user to select one
  if len(finv) == 1:
    finv = finv[0]

  else:
    print "Please choose expression for the inverse X = f^{-1}(Y)\n"
    for opt, expr in enumerate(finv):
      print "\nOption %d:\n----------\n" % opt
      sp.pprint(expr)
    choice = int(raw_input("Your Choice: "))
    finv = finv[choice]
    os.system('clear')

  return finv


# Function to output
def ItosLemma(f, drift, diffusion, finv=None, replaceX=True, pretty_print=False):
  """
    This function applies Ito's lemma given a twice-differentiable,
    invertible function (f) and drift and diffusion coefficients for a
    stochastic process (X_t) satisfying SDE

      dX = drift(X,t) dt + diffusion(X,t) dW

    This routine then prints out a new the SDE for Y=f(X)

      dY = drift_new(X,t) dt + diffusion_new(X,t) dW_t

    If replaceX=True, this will replace instances of X with f^{-1}(Y)

    NOTE: Inputs f, drift, diffusion _all_ expected to be sympy fcns.
    Derivatives are taken with respect to sympy symbols/variables X and
    t in those functions.
  """

  # Define symbols for display later
  t, Y, dt, dW, dY, dX, Y_t = sp.symbols('t Y dt dW dY dX Y_t')

  # Define differentiation functions
  DX, Dt = map(lambda var: (lambda fcn: sp.diff(fcn, var)), [X, t])

  # Compute drift_new and diffusion_new according to Ito's Lemma
  drift_new     = Dt(f) + DX(f)*b + 0.5*DX(DX(f))*(s**2)
  diffusion_new = DX(f)*s

  # If finv not given, invert Y=f(X) to get Y = f^{-1}(X)
  if finv == None:
    finv = invertInteractive(f, X, Y)

  # Substitute out X with f^{-1}(Y)
  if replaceX:
    drift_new, diffusion_new = map(lambda fcn: fcn.subs(X,finv).simplify(),
                                   [drift_new, diffusion_new])

  # Print the results to screen
  if pretty_print:
    sp.init_printing()
  print "\nOriginal SDE\n------------\n"
  sp.pprint(sp.Eq(dX, (drift)*dt + (diffusion)*dW))
  print "\n\nNew SDE for Y = f(X)\n--------------------\n"
  sp.pprint(sp.Eq(dY,(drift_new)*dt + (diffusion_new)*dW))
  print "\nwhere"
  sp.pprint(sp.Eq(Y,f))
  print("\n\n")


# Call function
ItosLemma(f, b, s)
