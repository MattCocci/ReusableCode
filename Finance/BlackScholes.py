##################################################
## Description ###################################
##################################################

# Simulates the price path of an asset S according
#   to the Black-Scholes SDE

# This script sets up two objects
#   1. A asset process object for S
#   2. A derivative object which defines the 
#	payoff


##################################################

import numpy as np
from matplotlib.pyplot as plt
from scipy import random.standard_normal


# Set up an object for the asset price process
class BlackScholesModel:	
    """ sigma	Volatility
	r	Risk-free rate
	d	Dividend yield
	S0	Initial asset price
    """
    
    # Store attributes
    def __init__(self, S0, sigma, r, d=0):
	self.S0 = float(S0)
	self.sigma = float(sigma)
	self.r = float(r)
	self.d = float(d)

    # Method to generate a simulation T into
    #	the future, N times
    def sim(self, T, N=1):
	
	# Draw N random normals
	W = standard_normal(N)

	# Compute the terminal value of the asset
	self.sim_path = self.S0 * np.exp( 
		T*(self.r - self.d - (1/2.) * self.sigma**2) 
		+ self.sigma * W)

    def plot(self):
	print 1


