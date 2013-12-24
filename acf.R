################################################################
## ACF #########################################################
################################################################

## Description ##
    # Assumes we have some data "y"
    # Computes the autocorrelation coefficient
    # Uses unbiased estimators for var and covar


autocorr.coef <- function(y, lags=1){

    # Ensures the process is demeaned
    if(mean(y) != 0){
	constant <- mean(y)
	y <- y - constant
    }



}





