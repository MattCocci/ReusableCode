
## Code for applying the vec operator to a matrix

vec <- function(A){
    return(cbind(as.vector(A)))
}

# Unstacks the vec operator into a matrix
#   Assumes vecA is a column matrix
unvec <- function(vecA, rows){

    # Rows is the number of rows in the matrix
    if(nrow(vecA) %% rows != 0){
	return("Length of vecA not evenly divisible
	       by given rows")
    } else {
	A <- matrix(vecA, rows, nrow(vecA)/rows)
	return(A)
    }
}

