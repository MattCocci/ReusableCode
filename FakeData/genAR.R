################################################################
## genAR #######################################################
################################################################

## Description ##
    # Generate AR(1) process which evolves according to 
    #	    y_t = rho * y_{t-1} + s*z
    #	then adds back in the mean, mu
    # rho: corr coeff
    # sd.e: standard deviation of the errors
    # n: length of series
    # mu: mean
    


# Set process parameters
rho  <- 0.8
sd.e <- 1 
n    <- 100
mu   <- 0


# Get the building blocks of the model
e    <- rnorm(n)
y    <- rep(0, n+1)


for(i in 1:n){
    y[i+1] = rho*y[i] + sd.e*e[i]
}


# generate series and lagged series
y.lag <- y[1:n] + mu
y <- y[2:(n+1)] + mu


