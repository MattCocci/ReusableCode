################################################################
## Regression ##################################################
################################################################

# This module creates all the functions you need to implement
#   regression in R because even though it's there already,
#   using preloaded functions just encourages doing
#   without understanding

# We use several results from linear algebra to speed up
#   computation

################################################################

# Regression of X on y where X and y are matrices
# NOTE: X does NOT include a constant already
# Dummify.cols indicates the columns with categorical
#   variables that we want to replace with dummy
#   variables; you set dummify.cols = c(x, y, ...)
#   when you run reg
reg <- function(X, y, cons=1, dummify.cols=NULL, homosked=0){

    ## Remove rows with an NA ##
    ind <- !(is.na(y) & rowSums(is.na(X) %*% diag(n)))
    X <- X[ind, ]
    y <- y[ind, ]


    ## Dummify categorical variables ##
    n.col <- ncol(X)
    if(dummify.cols != NULL){

	to.dum <- X[ , dummify.cols]
	X <- X[ , -dummify.cols]

	# Loop through columns to dummify
	for(col in 1:ncol(to.dum)){

	    # Extract column and column name
	    to.dum.col <- to.dum[ ,col]
	    name.stub <- colnames(to.dum)[col]

	    col.vals <- sort(unique(to.dum.col))

	    j = 1
	    for(val in col.vals){
		X <- cbind(X, to.dum.col==val)
		colnames(X)[n.col+j] <- paste(name.stub, j, sep="")
	    }
	}
    }


    # Observations
    n <- nrow(X) 
    k <- ncol(X)

    # Add a constant column if specified
    if(cons == 1){  	
	cbind(rep(1, n), X)
    }

    # Get coefficient estimate
    beta.hat <- solve(t(X) %*% X) %*% t(X) %*% y

    # Get the projection matrix and the orthogonoal 
    #	projection matrix and M.star
    # They allow for quick calculation of residuals and such
    P <- X %*% solve(t(X) %*% X) %*% t(X)
    M <- diag(n) - P
    M.star <- solve(diag(n) - (diag(n) * P))

    # Leverage values for all observations
    leverage <- (diag(n) * P) %*% matrix(rep(1, n), n, 1)

    # Fitted values (values fit using the entire sample
    #	of n observations)
    y.hat <- P %*% y

    # Residuals
    residuals <- M %*% y

    # Standardized residuals
    std.residuals <- (M.star^(1/2)) %*% residuals

    # Get the leave-one-out prediction errors 
    pred.error <- M.star %*% residuals

    ## Compute covariance matrix estimators for Beta##
    if(homosked==1){ # If forced homoskedastic errors
	s2 <- (1/(n-k)) * (t(residuals) %*% residuals)
	V <- (solve((1/n) * (t(X) %*% X)) * s2)/n

    } else { #DEFAULT CASE:
    # Compute ROBUST covariance matrix estimators

    # Eicker-White covariance matrix estimator
    D.hat <- diag(residuals^2)
    V.white <- solve((1/n) * (t(X) %*% X)) %*% ((1/n) * (t(X) %*% D.hat %*% X)) %*% solve((1/n) * (t(X) %*% X))
    V.white <- (n/(n-k))*V.white 
	# Rescaling so that it's unbiased (it's usually
	# biased towards zero)
    
    D.bar <- diag(std.residuals^2) 
    V.bar <- solve((1/n) * (t(X) %*% X)) %*% ((1/n) * (t(X) %*% D.bar %*% X)) %*% solve((1/n) * (t(X) %*% X))
	
    D.tilde <- diag(pred.error^2)
    V.tilde <- solve((1/n) * (t(X) %*% X)) %*% ((1/n) * (t(X) %*% D.tilde %*% X)) %*% solve((1/n) * (t(X) %*% X))

    # The standard output will be the zero-bias
    #	Horn, Horn, and Duncan cov matrix estimator
    V <- V.bar/n
    }

    # Get standard errors
    std.errors <- sqrt((V * diag(k)) %*% matrix(1, k, k))

    # Get the coefficient change from full sample estimate
    #	to leave-one-out estimate of beta
    coeff.chg <- solve(t(X) %*% X) %*% t(X) %*% (diag(pred.error))
    pred.chg <- (diag(n) * P) %*% pred.error
    

    # R squared
    R.squared <- 1 - (sum(residuals^2) / sum((y - mean(y))^2))

    # Adjusted R squared
    R.Squared.adj <- 1 - ((n-1)*(t(residuals) %*% residuals)) /
			((n-k)*(t(y-mean(y)) %*% (y-mean(y))))

    # R squared from prediction errors
    R.Squared.tilde <- 1 - (t(pred.error) %*% pred.error) /
			    (t(y-mean(y)) %*% (y-mean(y)))


    # Return the estimates
    toRet <- list("Coefficients" = beta.hat, 
		  "Standard Errors" = std.errors,
		  "Cov Matrix Estimator" = V, 
		  "Fitted Values" = y.hat,
		  "Residuals" = residuals,
		  "Standardized Residuals" = std.residuals,
		  "Predicted Error" = pred.error,
		  "Leverage Values" = leverage,
		  "Coefficient Change bc of Obs" = coeff.chg,
		  "Prediction Change bc of Obs" = pred.chg,
		  "R-Squared Tilde (BEST)" = R.squared.tilde,
		  "Adj. R-Squared" = R.squared.adj,
		  "R-Squared" = R.squared)
    return(toRet)
}

