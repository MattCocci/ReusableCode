################################################################
## VARAnalysis #################################################
################################################################

## Description ##
    # Example of VAR analysis
    # 

## IMPORT OTHER FUNCTIONS ######################################

library("MASS")
source("../kronecker.R")
source("../vec.R")

## IMPORT & FIX DATA ###########################################

data <- read.table("data.txt", header=TRUE)
dates <- seq(from=1959.0, to=2013.25, by=0.25)

# Get number of variables, subtracting one for dates
n <- ncol(data) - 1 # number of components/variables

# Log transformations
to.log <- c("rgdp", "wages", "m1")
for(i in 1:length(to.log)){
    data[to.log[i]] <- log(data[to.log[i]])
}

## IMPORT & FIX DATA ###########################################

# Y goes from start.obs to T+start.obs-1
# X goes from time start.obs-1 to T+start.obs-2
# Presample is everything before start obs
p <- 4 # lags
bigT <- 150  # size, T
start.obs <- 20 # Index where the response variable, Y, starts

# Build matrices to use in estimation
X <- matrix(0, bigT, n*p+1)
Y <- matrix(0, bigT, n)
for(lil.t in start.obs:(start.obs+bigT-1)){

    lil.x <- c(as.vector(t(data[(lil.t-1):(lil.t-p), 2:(n+1)])),1)
    X[lil.t-start.obs+1,] <- lil.x
    Y[lil.t-start.obs+1,] <- as.matrix(data[lil.t, 2:(n+1)])

}

# Get OLS estimator
phi.ols <- solve(t(X) %*% X) %*% t(X) %*% Y


## OLS FORECASTS ###########################################


    # Forecast from bigT + p + 1 onwards
    Y.forecast <- cbind(data$date, matrix(0, nrow(data), n))
    for(lil.t in 1:nrow(data)){

	# If before first forecast at time T+p+1, take historical
	if(lil.t <= bigT + start.obs - 1){  

	    Y.forecast[lil.t, 2:(n+1)] <- 
		as.matrix(data[lil.t, 2:(n+1)])

	} else {

	    # Else, forecast from last three lags
	    lil.x <- c(as.vector(
		 t(Y.forecast[(lil.t-1):(lil.t-p), 2:(n+1)])
		), 1)
	    Y.forecast[lil.t, 2:(n+1)] <- t(phi.ols) %*% lil.x
	}
    }

	
    # Plot the data
    titles = c("Log Real GDP", "Unemployment Rate", "Price Level",
	       "Log Wages", "M1")

    plot(dates, Y.forecast[,2], type="l", col="red", 
	 main=titles[1], xlab="", ylab=titles[1], lwd=2)
    lines(dates, data[,2], type="l", lwd=2)

    x11()

    par(mfrow = c(2, 2))
    for(i in 1:4){
	plot(dates, Y.forecast[,i+2], type="l", col="red", 
	     main=titles[i+1], xlab="", ylab=titles[i+1], lwd=2)
	lines(dates, data[,i+2], type="l", main=titles[i+1], 
	     xlab="", ylab=titles[i+1], lwd=2)
    }


    
## BAYESIAN FORECASTS ########################################





#
