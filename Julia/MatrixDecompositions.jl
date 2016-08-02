""" Slow Gaussian elimination of A to upper triangular form """
function gaussian_elimination_slow!(A::Array)
  n = size(A)[1]

  # col indexes column you are zeroing out below entry A[col,col]
  for col = 1:n-1
    # Current row below entry A[col,col] you are working on zeroing out
    for row = col+1:n
      multiplier    =  A[row,col]/A[col,col]
      A[row,col:n] += -multiplier*A[col,col:n]
      A[row,col]    = 0 # To truly zero out and avoid rounding error
    end
  end
  return A
end


"""
Given any matrix A with n rows, return the elimination matrices
E and M[:,:,1],...,M[:,:,n-1] where the following product returns an
upper triangular matrix

  E * A = (M[:,:,n-1] * ... * M[:,:,1]) * A

In the for-loop, we are just doing gaussian elimination steps at each
iteration, and storing each elimination step in a matrix. With the
resulting matrices, the product M[:,:,1]*A is equivalent to doing the
first Gaussian elimination step, M[:,:,2]*M[:,:,1]*A the first two
steps, etc.

Finally, E = M[:,:,n-1] * ... * M[:,:,1] is then just a stringing
together of all the elimination matrices into a single elimination
matrix E such that E*A is upper triangular.

Note: There is an option to do "pivoting," or swapping of rows in A so
that, in the process of Gaussian elimination, the multipliers aren't
huge. Particularly, if we're working on column "col", we find the
largest element A[i,col] where i >= col (i.e. below the diagonal), and
then swap the row A[col,:] with A[i,:]. That way, we avoid big
multipliers.
"""
function elimination_matrices(A_, pivot=false)
  n = size(A_)[1]    # Number of rows
  A = deepcopy(A_)   # Make a copy since we will be altering matrix
  R = zeros(n,n,n-1) # Will hold reordering matrices in event of pivot
  M = zeros(n,n,n-1) # Will hold elimination matrices
  E = eye(n)         # Product of all elimination and reordering matrices
  P = eye(Int32,n)   # Product of all reordering matrices

  # On any given iteration, we're zeroing out under A[col,col]
  # via Gaussian elimination
  for col = 1:n-1

    # Initialize reordering and elimination matrices
    R[:,:,col] = eye(n)
    M[:,:,col] = eye(n)

    # Before we do Gaussian elimination, check whether to pivot
    # (i.e. whether to swap rows so the multipliers aren't huge)
    if pivot
      # Get index of largest element in column "col" at/below diagonal
      ind_big = findmax( abs(A[col:n, col]) )[2][1]
      ind_big = ind_big + (col-1)
        #^Adjusts to proper index since findmax over col:n, not 1:n

      # Set reordering matrix by swapping rows
      R[ [col,ind_big], :, col] = R[ [ind_big,col], :, col]

      # Swap rows in A and p
      A = R[:,:,col]*A
    end

    # Do Gaussian elimination
    multipliers       = A[col+1:n,col] / A[col,col]
    A[col+1:n,col:n] += -multipliers * A[col,col:n]
    A[col+1:n,col]    = 0 # To truly zero out and avoid rounding error

    # Store gaussian elimination step in an elimination matrix
    M[col+1:n,col,col] = -multipliers

    # Update cumulative ordering & "reordering and elimination" matrices
    P = R[:,:,col]*P
    E = M[:,:,col]*R[:,:,col]*E
  end

  # Return elimination matrices and permutation vector
  return R, M, E, P
end


"""
Do Gaussian elimination on matrix A.

Note, A can be non-square. This allows for the case where there are more
columns than rows, which might happen if we want to solve system Ax=b
and passed [A b] as the argument to this function to row-reduce with b
included (then break up later).

There is an option to use pivoting in Gaussian elimination. In that
case, A[p,:] gives the proper row-ordering, so that E*A[p,:] delivers an
upper-triangular matrix.
"""
function gaussian_elimination(A::Array, pivot=false)
  R, M, E, P = elimination_matrices(A, pivot)
  return E*A, P
end



"""
Given square matrix A with n rows, return
- Lower triangular L
- Upper triangular U
- Reordering matrix
such that P*A = L*U

First, we get the reordering and elimination matrices, which satisfy

  E = M[:,:,n-1] * R[:,:,n-1] * ... * M[:,:,1] * R[:,:,1]
  P = R[:,:,n-1] * ... * R[:,:,1]

Next, note that, provided inv(E) exists (it does), clearly

  P*A = P*inv(E)*E*A

Since we know that E*A is upper triangular by construction, we will take

  U := E*A

Next, provided P*inv(E) is lower triangular, we can take

  L := P*inv(E) = P * inv(R[:,:,1])   * inv(M[:,:,1])   * ... *
                      inv(R[:,:,n-1]) * inv(M[:,:,n-1])

But note that both the R and M matrices have special forms. Since the R
matrices are just identity matrices with rows swapped, it is the case
that R[:,:,i] = inv(R[:,:,i]) for all i. Therefore,

  L := P*inv(E) = P * R[:,:,1]   * inv(M[:,:,1])   * ... *
                      R[:,:,n-1] * inv(M[:,:,n-1])

Next, each M[:,:,i] is lower triangular, with only column i having
nonzero entries below the diagonal. We can form the inverse simply by
negating all elements below the diagonal (i.e. all elements in column i
below the diagonal).  Finally, writing out P and simplifying, it is the
case that L indeed ends up being lower triangular.
"""
function LU_decomposition(A, pivot=true)
  R, M, E, P = elimination_matrices(A, pivot)

  n = size(A)[1]
  U = E*A
  L = P
  for i = 1:n-1
    # Start with M[:,:,i], form inverse by negating elements in col i,
    # below row i
    Minv = M[:,:,i]
    Minv[i+1:n,i] = -Minv[i+1:n,i]

    # L = D * R[:,:,1] * inv(M[:,:,1]) * ... * R[:,:,n-1] * inv(M[:,:,n-1])
    L = L*R[:,:,i]*Minv
  end
  return L, U, P
end



""" Forward-solve system Lx = b where L lower triangular """
function solve_forward(L,b)
  x = deepcopy(b)/L[1,1]
  for i = 2:n
    x[i] = ( b[i] - dot(vec(L[i,1:i-1]), x[1:i-1]) ) / L[i,i]
  end
  return x
end


""" Back-solve system Ux = b where U upper triangular """
function solve_backward(U,b)
  n = length(b)
  x = deepcopy(b)/U[n,n]
  for i = (n-1):-1:1
    x[i] = ( b[i] - dot(vec(U[i,i+1:n]), x[i+1:n]) ) / U[i,i]
  end
  return x
end



"""
Solve system Ax = b like A \ b

We do this as Matlab does. First, the LU decomposition gives us P, L, U
such that PA = LU. Then, we can write the system Ax = b as

  PAx = LUx = Pb

To solve, define c := Ux. Then solve (i) forward for c, then (ii)
backward for x

  (i) Lc = Pb    (ii) Ux = c

This is tantamount to computing x = inv(U)*inv(L)*P*b. It's just that,
after the LU decomposition, we're solving simple triangular systems.
"""
function solve_system(A,b)

  # Do an LU decomposition
  L, U, P = LU_decomposition(A)

  c = solve_forward(L,P*b)
  x = solve_backward(U,c)
  return x
end



"""
Cholesky Factorization of Symmetric Positive Definition Matrix A = CC'
where C is lower triangular.

We solve recursively

  A = [a1  a2'] = [c1 0] [c1 c2']  = CC'
      [a2  B  ]   [c2 D] [0  D' ]

where a and l are scalars, b and ll are vectors, and A and L matrices.
This gives system of equations

  a1 = (c1)^2         ->    c1   = sqrt(a1)
  a2 = c1*c2          ->    c2   = a2/c1
  B  = c2*c2' + D*D'  ->    D*D' = B - c2*c2'

"""
function cholesky_decomposition(A)
  # Partition A matrix
  a1, a2, B  = A[1,1], A[2:end,1], A[2:end,2:end]

  # Cholesky computations
  c1 = sqrt(a1)
  c2 = a2/c1

  # If the A matrix that came in was 2x2, you are at the last step.
  # Can compute D and conclude
  if size(A)[1] == 2
    D = sqrt(B - c2*c2')
    return [c1 0; c2 D]

  # Otherwise, recurse
  else
    D = cholesky_decomposition(B - c2*c2')
    return [c1 zeros(1,length(c2)); c2 D]
  end
end


"""
Create Householder Matrix from vector u:

  H = I - 2*u*u' / ||u||^2
"""
HouseholderMatrix(u) = ( eye(length(u)) - 2*(u*u')/(norm(u)^2) )



"""
Multiply A by Householder matrix associated with vector u

Given u, directly compute

  A - 2*u*u'*A / ||u||^2

which is equivalent to H*A, where H is the Householder matrix
associated with vector u. This is just computationally quicker than
forming H and computing H*A.
"""
HouseholderMultiply(u,A) = A - 2*(u*u'*A)/(norm(u)^2)



"""
QR Decomposition of a Matrix

Given A, we compute orthogonal Q and upper triangular R such that

  A = Q*R

We construct R by a sequence of orthogonal Householder updates to A.
In particular, start with matrix (m x n) A. We construct the Householder
matrix from the vector u, constructed as follows

  a     = A[:,1]      # First column of A
  alpha = +- norm(a)  # With the sign the opposite of element A[1,1]
  u     = a - alpha*[1; 0; ...; 0]

From u, construct H. It can be shown that H*A will have all zeros under
the first element of column 1. We can repeat this process, constructing
a sequence of Householder matrices for lower-right blocks that haven't
yet been put into upper triangular form yet.

Then, if H_1,...,H_p for p=min(m,n) is the sequence of Householder
matrices used to get A into upper triangular form, we have

  R = H_p * ... * H_1 * A

hence

  Q = inv(H_p * ... * H_1) = inv(H_1) * ... * inv(H_p)

Since Householder matrices are orthogonal and symmetric, we have
inv(H_i) = H_i' = H_i so

  Q = H_1 * ... * H_p

"""
function QR_decomposition(A)
  m,n  = size(A)
  R    = deepcopy(A)
  Q    = eye(m)
  I    = eye(m)
  stop = min(m,n)

  # Loop over columns
  for i = 1:stop
    # Grab column i, zeroing out everything above R[i,i].
    a = [zeros(i-1,1); R[i:m,i]]

    # Construct vector and associated Householder matrix, which will
    # reduce R[i:m,i:n] and preserve reductions in R[1:i-1,1:i-1]
    alpha = -sign(a[i])*norm(a)
    u     = a - alpha*I[:,i]
    H     = HouseholderMatrix(u)

    # Householder update to R
    R = HouseholderMultiply(u, R)

    # Update Q
    Q = Q*H

    # zero out below diagonal for good measure
    R[i+1:m,i] = 0
  end

  return Q, R
end



"""
Given vector y and matrix X, return least squares solution to

  y = X*b + e

We think of this as an overdetermined system Xb ~= y, compute the QR
decomposition of X, and define

  Q_1 = Q[1:m,1:n]
  R_1 = R[1:n,1:n]

This is sufficient since the regression coefficients are simply the
solutions to the upper-triangular least squares problem

  R_1*x = Q_1'*b

"""
function QR_leastsquares(y,X)
  m, n = size(X)
  Q, R = QR_decomposition(X)
  Q_1  = Q[1:m,1:n]
  R_1  = R[1:n,1:n]
  return solve_backward(R_1, Q_1'*y)
end
