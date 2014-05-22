
# Bisection for root-finding

def bisect(f, lower, upper, tol=10e-5):

    middle = 0.5 * (upper + lower)
    print "Middle " + str(middle)

    if upper - lower >= tol:
	if f(middle) > 0:
	    bisect(f, lower, middle, tol)
	else:
	    bisect(f, middle, upper, tol) 
    else:
	return middle
