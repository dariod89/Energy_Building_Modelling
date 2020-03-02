# This function computes the correlation between NxM pairs of k-dimensional inputs. 
# Correlation lengths must be provided in vector d. The input 'string' specifies which correlation function to use.
#
# INPUTS:
# X:      Nxk matrix of N k-dim vectors
# Y:      Mxk matrix of N k-dim vectors
# d:      k-dim vector of positive correlation lengths
# string: a string specifying which correlation function to use.
#         One of: 'exp2', ... (others will be added later)
#
# OUTPUT:
# A: NxM matrix, A[i,j] is the correlation between vectors X[i,] and Y[j,]
#
# The function uses multidimensional product (between tensors of order 3) to increase speed.
# The choice necessarily makes the code less straightforward to interpret. However, what it does is really simple:
#
# Starting from X, Y and d as above, the element (i,j) of the matrix A is 
#
#                  A[i,j]=exp(-S), 
#
# where S is the sum of the squares of the k components of the vector z = (X[i,]-Y[,j])./d  
# ("./" denotes division component by component). Ie, S = sum of z[h]^2  for h=1, ...k.
#
# A slower but straightforward-to-interpret version of the same function, implemented via two nested for loops,
# is provided after the following function.

corr.fun <- function(X, Y, d, string){

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  d <- as.numeric(d)
  
  N <- dim(X)[1]
  M <- dim(Y)[1]
  k <- dim(X)[2]
  
  # Create X3 and Y3 NxMxk arrays, by stacking up original X and Y 
  X3 <- array(X, dim=c(dim(X),M) )  # NxkxM
  X3 <- aperm(X3, c(1,3,2))         # NxMxk
  Y3 <- array(Y, dim=c(dim(Y),N))   # MxkxN
  Y3 <- aperm(Y3, c(3,1,2))         # NxMxk  
  
  Z <- X3-Y3  # NxMxk
  for (i in 1:k)
    Z[,,i] <- Z[,,i]/d[i] # rescale dimensions by corr length
  
  if (identical(string,'exp2')){
    A <- rowSums(Z^2, dims=2) # A[i,j] = sum of the elements Z[i,j,:]^2 
    A <- exp(-A/2)
  }
  return(A)
  
}


# Same function as above, implemented via nested for loops: straightforward to interpret, but slower.
corr.fun2 <- function(X, Y, d, string){

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  N <- dim(X)[1]
  M <- dim(Y)[1]
  k <- dim(X)[2]
  
  D <- diag(1/d)
  A <- matrix(0, N,M)
  for (i in 1:N)
    for (j in 1:M){
      x = D%*%(X[i,]-Y[j,])
      x = sum(x^2)
      A[i,j] <- exp(- x/2)
    }
  return(A) 
}


# This function carries out leave-one-out cross validation over an emulator. 
# In turn, it leaves out one observation, builds the emulator without it, and returns the emulator 
# predicted distribution (mean and variance) for the output corresponding to the  left-out point.

Cross_Val <- function(Design.par.full, Design.regr.full, y.full, beta, Cov.beta, d, nu, sigma2, string){
  N <- dim(Design.par.full)[1]
  p <- dim(Design.par.full)[2]
  q <- dim(Design.regr.full)[2]
  
  ind.full <- 1:N
  
  M <- array(dim=c(N,1))
  V <- array(dim=c(N,1))
  
  for (i in 1:N){
    ind <- ind.full[-i]
    Design.par  <- Design.par.full[ind, , drop=FALSE]
    Design.regr <- Design.regr.full[ind, , drop=FALSE]
    test.par    <- Design.par.full[i, , drop=FALSE]
    test.regr   <- Design.regr.full[i, , drop=FALSE]
    y <- y.full[ind]
    
    res <- emul(test.par, test.regr, Design.par, Design.regr, y, beta, Cov.beta, d, nu, sigma2, string)
    M[i] <- res$Mean
    V[i] <- res$Variance
  }
  
  return( list("Mean"=M, "Variance"=V) )
}







