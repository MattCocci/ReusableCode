import sympy as sp
import os

# Initialize symbols
X, mu, sigma = sp.symbols('X mu sigma')

# Geometric Brownian Motion
b = mu*X
s = sigma*X
f = sp.ln(X)
outfile = 'GBM.tex'

# Some Other Example
#b = 0
#s = 1
#f = X**2
#outfile = False


# Function to output
def ItosLemma(f, b, s, outfile=False, pretty_print=False):
  """
    This function applies Ito's lemma given a twice-differentiable,
    invertible function (f) and drift and diffusion coefficients for a
    stochastic process (X_t) satisfying SDE

      dX = b(X,t) dt + s(X,t) dW

    Inputs are _all_  expected to be sympy functions with sympy
    symbols/variables X (capitalized) and/or t in them.

    This routine then prints out a new the SDE

      dY = b*(Y,t) dt + s*(Y,t) dW_t

    where Y = f(X)
  """

  # Define symbols for display later
  t, Y, dt, dW, dY, dX, Y_t = sp.symbols('t Y dt dW dY dX Y_t')

  # Compute b* and s* according to Ito's Lemma
  bstar = sp.diff(f,t) + sp.diff(f,X)*b + 0.5*sp.diff(f,X,X)*(s**2)
  sstar = sp.diff(f,X)*s

  # Make bstar and sstar a function of y=f(x) and t, not _x_ and t
  finv = sp.solve(Y-f,X)  # Calculate inverse of, i.e. get x = f^{-1}(y)
  if len(finv) > 1:
    print "Please choose expression for inverse X = f^{-1}(Y)\n"
    for opt, expr in enumerate(finv):
      print "\nOption %d:\n----------\n" % opt
      sp.pprint(expr)
    choice = int(raw_input("Your Choice: "))
    finv = finv[choice]
    os.system('clear')
  else:
    finv = finv[0]

  bstar = bstar.subs(X,finv).simplify() # Substitute x = f^{-1}(y)
  sstar = sstar.subs(X,finv).simplify() # Substitute x = f^{-1}(y)

  # Print the results to screen
  if pretty_print:
    sp.init_printing()
  print "\nOriginal SDE\n------------\n"
  sp.pprint(sp.Eq(dX,(b)*dt + (s)*dW))
  print "\n\nNew SDE for Y = f(X)\n--------------------\n"
  sp.pprint(sp.Eq(dY,(bstar)*dt + (sstar)*dW))
  print "\nwhere"
  sp.pprint(sp.Eq(Y,f))
  print("\n\n")

  # Build up latex string
  if outfile:
    tex = 'dY_t = ('
    tex += sp.latex(bstar.subs(Y,Y_t)) + ')\\; dt + ('
    tex += sp.latex(sstar.subs(Y,Y_t)) + ')\\; dW_t'
    f = open(outfile, 'w')
    f.write('\\begin{equation}\n'+ tex + '\n\\end{equation}')
    f.close()


# Call function
ItosLemma(f, b, s, outfile=outfile, pretty_print=True)
