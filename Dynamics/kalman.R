################################################################
## Kalman Filter ###############################################
################################################################

# This module creates all the functions you need to implement
#   a Kalman Filter

# This is based off of kalman.py, the python module crated by
#   Stachurski and Sargent for their Quant-Econ site


# Implements the Kalman filter for the state space model
#
#    x_{t+1} = A x_t + w_{t+1}
#    y_t = G x_t + v_t.
#
# Here x_t is the hidden state and y_t is the measurement.  
# The shocks {w_t} and {v_t} are iid zero mean Gaussians with 
# covariance matrices Q and R respectively.

################################################################


# Set up your prior, where x is a column matrix; Sigma a matrix
x0 <- #matrix(vector) 
Sigma0 <- #matrix() 

# Set up the law of motion




