#' Compute the solution of weighted LS + fused SCAD penalty. 
#'
#' @param y The response vector.
#' @param X The design matrix.
#' @param w The vector of weights in the least squares term.
#' @param D The fusion matrix.
#' @param lambda The regularization parameter.
#' @param a The hyper-parameter in SCAD.
#' @param beta0 The initial vector of coefficients.
#' @param maxiters The maximum number of iterations for the algorithm
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
fused_SCAD <- function(y, X, w=rep(1, length(y)), D, lambda, a=3.7, beta0=rep(0, ncol(X)), maxiters=1e2, tol=1e-5) {
  # the key is that we transform it to a fused lasso problem

  for (k in 1:maxiters) {
    # for each iteration, D changes
    Dbeta <- D%*%beta0
    L <- Diagonal(n=nrow(D), x=scad_prime(Dbeta, lambda=lambda, a=a))
    D1 <- L%*%D
    
    # at each iteration, solve a fused lasso problem
    beta <- fused_Lasso(y, X, w, D=D1, lambda=1)$beta #Note here lambda=1 for the fused lasso problem
    
    #check convergence
    if ( norm(beta-beta0,'2') < tol*(1 + norm(beta0,'2')) ) break
   
    # update the estimate
    beta0 <- beta
  }
  
  return(list(beta=beta, iters=k))
}






#' Compute the derivative of SCAD penalty.
#'
#' @param beta The vector of coefficients.
#' @param lambda The regularization parameter.
#' @param a The hyper-parameter in SCAD.
#' @export
scad_prime <- function(beta, lambda, a=3.7){
  prime <- rep(0, length(beta))
  for (j in 1:length(beta)) {
    prime[j] <- ifelse(abs(beta[j]) <= lambda, lambda, max(a*lambda-abs(beta[j]), 0)/(a-1))
  }
  
  return(prime)
}



