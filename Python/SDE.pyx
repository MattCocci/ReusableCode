#import sympy as sp
#import numpy as np
#import warnings

print "Hello Word"

#class Discretize:
  #def __init__(self, b, s, soln=None):

    ## Store symbolic functions and derivatives
    #for diffs in range(3):
      #prfx = diffs*'d'
      #setattr(self, prfx+'b_sym', b)
      #setattr(self, prfx+'s_sym', s)
      #b, s = map(lambda f: sp.diff(f,'X'), [b, s])

    #for diffs in range(3):
      #prfx = diffs*'d'
      #for atr in ['b', 's']:
        #setattr(self, prfx+atr, \
          #sp.lambdify(X, getattr(self, prfx+atr+'_sym'), 'numpy'))

    #if soln != None and not isinstance(soln, type(lambda:0)):
      #soln = sp.lambdify((x0, t, W), soln, 'numpy')
    #self.soln = soln

  #def EM(self, Xt, dt, dW, dZ=None):
    #return Xt + self.b(Xt)*dt + self.s(Xt)*dW

  #def Milstein(self, Xt, dt, dW, dZ=None):
    #return ( self.EulerStep(Xt, dt, dW)
              #+ 0.5*self.s(Xt)*self.ds(Xt)*(dW**2-dt) )

  #def Taylor15(self, Xt, dt, dW, dZ):
    #return (
      #self.MilsteinStep(Xt, dt, dW)
      #+ (0.5)*(self.b*self.db + 0.5*(self.s**2)*self.ddb)*(dt**2)
      #+ self.s * self.db * dZ
      #+ (self.b*self.ds + 0.5*(self.s**2)*self.dds)*(dt*dW-dZ)
      #+ 0.5*((self.s**2)*self.dds + self.s*(self.ds**2))
        #* (1/3.)*((dW**2)-dt)*dW
      #)

  #def SimulatePaths(self, methods, x0, t0, T, dt, analytical=False):

    ## If the user wants an anlytical solution but didn't provide one at
    ## initialization or subsequently, error out
    #if analytical:
      #assert self.soln != None, \
        #'No analytical solution given at initialization; cannot compute.'
    #analytical = int(analytical)

    ## Initialize parameters and X matrix to hold simulated paths
    #P = np.arange(t0, T, dt) # Left endpoints of subintervals
    #nsteps   = len(P)        # Number of subintervals/steps we'll take
    #nmethods = len(methods)  # Number of methods
    #X = np.vstack( (x0*np.ones(nmethods+analytical),
                    #np.nan*np.ones((nsteps,nmethods+analytical))) )

    ## Generate random draws of just dWs (if EM and/or Milstein) or dWs
    ## and dZ if also doing Taylor 1.5
    #if 'Taylor15' in methods:
      #mean = np.zeros(2)
      #cov  = np.array([ [dt,     0.5*dt],
                        #[0.5*dt, (1/3.)*(dt**3)] ])
    #else:
      #mean = 0*np.ones(1)
      #cov  = dt*np.ones((1,1))

    ## Draws and functions to access them
    #draws = np.random.multivariate_normal(mean=mean, cov=cov, size=nsteps)
    #dW = lambda t: draws[t,0]
    #if 'Taylor15' in methods:
      #dZ = lambda t: draws[t,1]
    #else:
      #dZ = lambda t: None

    #for t in range(nsteps):
      #X[t+1,0:-(1+analytical)] = map(lambda m:
                      #getattr(self, methods[m])( X[t,m], dt, dW(t), dZ(t) ),
                      #range(len(methods)))
    #return X




