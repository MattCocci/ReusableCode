import sympy as sp
import numpy as np
import SDE_PurePython
import sys

X = sp.symbols("X")
b = -1*X
s = 1*X

D = SDE_PurePython.Discretize(b, s)
X = D.SimulatePaths(['EM', 'Milstein'], 1, 0, 10, float(sys.argv[1]))

