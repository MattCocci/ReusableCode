import sympy as sp
import warnings
import numpy as np
# Get special cpython/numpy stuff
cimport numpy as np

# Datatype for our arrays
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

X, x0, t, W = sp.symbols('X x0 t W')
class Discretize:
  def __init__(self, b, s, soln=None):

    # Store symbolic functions and derivatives
    for diffs in range(3):
      prfx = diffs*'d'
      setattr(self, prfx+'b_sym', b)
      setattr(self, prfx+'s_sym', s)
      b, s = map(lambda f: sp.diff(f,'X'), [b, s])

    for diffs in range(3):
      prfx = diffs*'d'
      for atr in ['b', 's']:
        setattr(self, prfx+atr, \
          sp.lambdify(X, getattr(self, prfx+atr+'_sym'), 'numpy'))

    if soln != None and not isinstance(soln, type(lambda:0)):
      soln = sp.lambdify((x0, t, W), soln, 'numpy')
    self.soln = soln

  def EM(self, float Xt, float dt, float dW, float dZ):
    return Xt + self.b(Xt)*dt + self.s(Xt)*dW

  def Milstein(self, float Xt, float dt, float dW, float dZ):
    return ( Xt + self.b(Xt)*dt + self.s(Xt)*dW
              + 0.5*self.s(Xt)*self.ds(Xt)*(dW**2-dt) )

  def Taylor15(self, float Xt, float dt, float dW, float dZ):
    return (
      self.MilsteinStep(Xt, dt, dW)
      + (0.5)*(self.b*self.db + 0.5*(self.s**2)*self.ddb)*(dt**2)
      + self.s * self.db * dZ
      + (self.b*self.ds + 0.5*(self.s**2)*self.dds)*(dt*dW-dZ)
      + 0.5*((self.s**2)*self.dds + self.s*(self.ds**2))
        * (1/3.)*((dW**2)-dt)*dW
      )

  def SimulatePaths(self, methods, float x0, float t0, float T, float dt, int analytical=0):

    # If the user wants an anlytical solution but didn't provide one at
    # initialization or subsequently, error out
    if analytical:
      assert self.soln != None, \
        'No analytical solution given at initialization; cannot compute.'

    # Initialize parameters and X matrix to hold simulated paths
    cdef np.ndarray P = np.arange(t0, T, dt, dtype=DTYPE) # Left endpoints of subintervals
    cdef int nsteps   = len(P)        # Number of subintervals/steps we'll take
    cdef int nmethods = len(methods)  # Number of methods
    cdef np.ndarray X = np.vstack( (x0*np.ones(nmethods, dtype=DTYPE),
                    np.nan*np.ones((nsteps,nmethods), dtype=DTYPE)) )

    # Determine if we'll need to draw another RV for Taylor 1.5 method
    cdef int dim = 1
    if 'Taylor15' in methods:
      dim = 2

    # Set up mean and cov vectors for draws
    cdef np.ndarray mean = np.zeros(dim, dtype=DTYPE)
    cdef np.ndarray cov  = np.zeros((dim,dim), dtype=DTYPE)
    cdef int zind = 0
    if dim > 1:
      cov  = np.array([ [dt,     0.5*dt],
                        [0.5*dt, (1/3.)*(dt**3)] ])
      zind = 0
    else:
      cov  = dt*np.ones((1,1))

    # Generate random draws of just dWs (if EM and/or Milstein) or dWs
    # and dZ if also doing Taylor 1.5
    cdef np.ndarray draws = np.random.multivariate_normal(mean=mean, cov=cov, size=nsteps)

    cdef int t, m
    for t in range(nsteps):
      X[t+1,:] = map(lambda m:
                      getattr(self, methods[m])( X[t,m], dt, draws[t,0], draws[t,zind]),
                      range(len(methods)))
    return X




