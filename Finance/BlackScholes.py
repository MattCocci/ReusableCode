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
import matplotlib.pyplot as plt
from scipy.stats import norm


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

        # Store the number of simulations
        self.sim_N = N
        
        # Draw N random normals
        W = norm.rvs(loc=0, scale=1, size=N)

        # Compute the terminal value of the asset
        self.sim_path = self.S0 * np.exp(float(T)*(self.r - self.d - (1/2.) * self.sigma**2) + self.sigma * W)

    # Method to plot a histogram of terminal values
    def plot(self):
        plt.hist(self.sim_path, normed=True, alpha=0.75)
        title = r'Distribution of Terminal Values, $S_0$ = ' + str(self.S0) + r', N = ' + str(self.sim_N)
        plt.title(str(title))
        plt.grid(True)
        plt.xlabel(r'Terminal Values, $S_T$')
        plt.ylabel('Probability')
        plt.show()

	

