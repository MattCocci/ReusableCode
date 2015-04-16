import sympy as sp
import numpy as np
import pyximport;
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True)
import SDE

X = sp.symbols("X")
b = -1*X
s = 1

D = SDE.Discretize(b, s)
X = D.SimulatePaths(['EM', 'Milstein'], 1, 0, 10, 0.000001)
