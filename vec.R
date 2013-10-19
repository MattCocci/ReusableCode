
## Code for applying the vec operator to a matrix

vec <- function(A){
    n <- ncol(A)

    toRet <- rep(0, 0)
    for(i in 1:n){  
	toRet <- c(toRet, A[,i]
    }
    return toRet
}

