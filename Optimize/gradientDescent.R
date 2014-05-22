

sigmoid <- function(theta, x){
    expon <- exp(theta*x)
    toRet <- expon/(1+expon)
    return(toRet)
}


gradientDescent <- function(X,y,theta, alpha, num.iters){
    
    # number of observations
    m <- length(y)
    J.history <- rep(0,num.iters)

    for(iter in 1:num.iters){

	theta = 

    }
    J.history[iter] <- computeCost(X,y,theta)
}



