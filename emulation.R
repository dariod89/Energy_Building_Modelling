# INPUTS
#
# test.par:    Nxp matrix/data.frame of raw input parameters at which to evaluate the emulator
# test.regr:   Nxq matrix/data.frame of regressors associated to parameters in 'test.par'
# Design.par:  nxp design matrix/data.frame of input parameters
# Design.regr: nxq design matrix/data.frame of regressors
# y:           n-dim vector of simulator responses
# beta:        q-dim vector of prior coefficients for regression surface
# Cov.beta:    qxq prior covariance matrix of beta
# d:           p-dim vector of correlation lengths on inputs
# nu:          variance of nugget process
# sigma2:      variance of the GP
# string:      specifies corr function to use. One of: 'exp2', 


emul <- function(test.par, test.regr, Design.par, Design.regr, y, beta, Cov.beta, d, nu, sigma2, string){
  T <- nrow(test.par)
  n <- nrow(Design.par)
  
  # PRIOR MEAN AND VARIANCE OF TEST INPUTS
  prior.mean <- test.regr %*% beta        # Tx1, prior mean at test.par
  a <- array(dim=c(T,1))                  # Tx1
  K <- Cov.beta %*% t(test.regr)          # qxT
  for (i in 1:T){
    vec <- test.regr[i,]
    a[i] <- vec %*% K[,i] # 1xq x qx1
  }
  prior.var <- a + sigma2 + nu            # Tx1, prior variance at test.par
  
  # PRIOR MEAN AND COVARIANCE OF DESIGN POINTS
  prior.D <- Design.regr %*% beta                          # n-dim vector
  a <- Design.regr %*% Cov.beta %*% t(Design.regr)         # nxn, prior covariance of regression part
  b <- sigma2*corr.fun(Design.par, Design.par, d, string)  # nxn, prior GP cov
  c <- nu*diag(n)
  A <- a + b + c                                           # nxn, prior covariance of design points
  
  # PRIOR COVARIANCE BETWEEN INPUTS AND DESIGN POINTS
  a <- test.regr %*% Cov.beta %*% t(Design.regr)         # Txn, regression part
  b <- sigma2*corr.fun(test.par, Design.par, d, string)  # Txn, prior GP covariance
  tx <- a + b                                            # Txn, vector t(x)
  
  # POSTERIOR MEAN AND COVARIANCE FOR INPUT PARAMETERS
  e <- y - prior.D                                # n-dim vector
  post.mean <- prior.mean + tx%*%solve(A,e)       # Tx1 vector
  a <- array(dim=c(T,1))                          # Tx1 vector
  K <- solve(A, t(tx))                            # nxT
  for (i in 1:T){
    vec <- tx[i,]
    a[i] <- vec %*% K[,i]
  }
  post.var  <- prior.var  - a   # NxN matrix
  
  return_list <- list("Mean"=post.mean, "Variance"=post.var)
  return(return_list)
}
