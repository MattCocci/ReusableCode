import sympy as sp
execfile("DiscretizeSDE.py")

alpha, X, sigma = sp.symbols("alpha X sigma")
b = -alpha*X
s = sigma

D = DiscretizeSDE(b, s)
X = D.SimulatePaths(['EM'], 1, 0, 10, 0.1)
