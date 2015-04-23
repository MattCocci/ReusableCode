# Setup file to compile SDE.pyx
#
# To compile, run this at the command line:
#
#   python setupSDE.py build_ext --inplace
#
from distutils.core import setup
from distutils.extension import Extension


# Whether to compile the C file or use Cython
USE_CYTHON = 0

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension('SDE', ['SDE'+ext])]

if USE_CYTHON:
  from Cython.Build import cythonize
  import numpy
  extensions = cythonize("SDE.pyx", include_path = [numpy.get_include()]) 

setup(
    ext_modules = extensions
)
