################################################################
## genLM #######################################################
################################################################

## Description ##
    # Generate Linear model (model with a linear CEF)
    #	    y_t = alpha + beta1 * x1 + beta2 * x2 + s.e * e
    #
    # (y_i, X_i)    Jointly normally distributed
    # alpha	    constant term
    # x1	    Normal continuous variable, centered at 0
    # x2	    Indicator variable; binomial and takes value
    #			of 1 with probability theta
    # theta	    The probability that x2 = 1
    # s.e	    Standard Error for normally distributed errors
    # n		    Size of the sample


################################################################
genLM <- function(alpha=5, beta1=2, beta2=10, theta=0.5, s.e=2, n=100){

    x1 <- rnorm(n)
    x2 <- rbinom(n, 1, theta)
    e  <- rnorm(n)
    y  <- alpha + beta1*x1 + beta2*x2 + s.e*e 

    y  <- matrix(y, n, 1) 
    X  <- cbind(x1, x2)
    colnames(y) <- "Income"
    colnames(X) <- c("Ability", "Female")

    toRet <- list("y" = y, 
		  "X" = X)
    return(toRet)

}
