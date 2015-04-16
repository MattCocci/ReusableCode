import sympy as sp
import numpy as np
import SDE_PurePython

X = sp.symbols("X")
b = -1*X
s = 1

D = SDE_PurePython.Discretize(b, s)
X = D.SimulatePaths(['EM', 'Milstein'], 1, 0, 10, 0.00001)

