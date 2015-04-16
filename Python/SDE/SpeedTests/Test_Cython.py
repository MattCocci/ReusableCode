import sympy as sp
import numpy as np
import sys
import pyximport;
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True)
import SDE

X = sp.symbols("X")
b = -1*X
s = 1*X

D = SDE.Discretize(b, s)
X = D.SimulatePaths(['EM','Milstein'], 1, 0, 10, float(sys.argv[1]))
