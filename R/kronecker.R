
# Calculates the Kronecker product
kronecker <- function(A, B){

    m <- nrow(A)
    n <- ncol(A)

    p <- nrow(B)
    q <- ncol(B)

    toRet <- matrix(0, m*p, n*q)
    for(i in 1:m){
	for(j in 1:n){
	    start.row <- 1+(i-1)*p
	    stop.row  <- i*p
	    start.col <- 1+(j-1)*q
	    stop.col  <- j*q
	    toRet[start.row:stop.row, start.col:stop.col] <- 
		A[i, j] * B
	}
    }
    return(toRet)

}

