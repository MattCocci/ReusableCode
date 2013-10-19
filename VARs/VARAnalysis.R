################################################################
## VARAnalysis #################################################
################################################################

## Description ##
    # Example of VAR analysis
    # 

## IMPORT & FIX DATA ###########################################

data <- read.table("data.txt", header=TRUE)
dates <- seq(from=1959.0, to=2013.25, by=0.25)

# Get number of variables, subtracting one for dates
n <- ncol(data) - 1 # number of components/variables

# Log transformations
to.log <- c("rgdp", "prices", "wages", "m1")
for(i in 1:length(to.log)){
    data[to.log[i]] <- log(data[to.log[i]])
}


## IMPORT & FIX DATA ###########################################

# Y goes from 1+p to T+p; X goes from time 1 to T+p-1
p <- 4 # lags
bigT <- 150  # size of same, T
last.obs <- bigT + p + 1

# Build matrices to use in estimation
X <- matrix(0, bigT, n*p+1)
Y <- matrix(0, bigT, n)
for(i in 1:bigT){

    lil.x <- c(as.vector(t(data[(i+p-1):i, 2:(n+1)])),1)
    X[i,] <- lil.x
    Y[i,] <- as.matrix(data[i+p, 2:(n+1)])

}


# Get OLS estimator
phi.ols <- solve(t(X) %*% X) %*% t(X) %*% Y

# Forecast from bigT + p + 1 onwards
Y.forecast <- cbind(data$date, matrix(0, nrow(data), n))
for(i in 1:nrow(data)){

    # If before first forecast at time T+p+1, take historical
    if(i <= bigT + p ){  

	Y.forecast[i, 2:(n+1)] <- as.matrix(data[i, 2:(n+1)])

    } else {

	# Else, forecast from last three lags
	lil.x <- c(as.vector(t(Y.forecast[(i-1):(i-p), 2:(n+1)])), 1)
	Y.forecast[i, 2:(n+1)] <- t(phi.ols) %*% lil.x
    }
}


## OLS FORECASTS ###########################################

    # Compare OLS forecasts to actual data
    #series <- 2 # rgdp
    #series <- 3 # unemp
    #series <- 4 # prices
    #series <- 5 # wages 
    series <- 6 # m1

    # Plot the forecast 
    plot(dates, Y.forecast[,series], type="l", col="red")

    # Plot the Actual data
    lines(dates, data[,series], col="black")
    
    






