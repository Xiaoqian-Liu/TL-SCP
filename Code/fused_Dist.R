#' Compute the solution of weighted LS + fusion-distance penalty.
#'
#' @param y The vector of response.
#' @param X The design matrix.
#' @param w The vector of weights in the least squares term.
#' @param D The fusion matrix.
#' @param k The regularization parameter k for the distance penalty (equal to the number of nonzeros in Dbeta). 
#' @param rho The initial value of rho, increases at each iteration to reach to a large value.
#' @param beta0 The initial vector of coefficients.
#' @param maxiters The maximum number of iterations.
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
fused_Dist <- function(y, X, w=rep(1, length(y)), D, k, rho=1,  beta0=rep(0, ncol(X)), maxiters=1e2, tol=1e-5){
  n <- length(w)
  W <- Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X
  beta <- beta0
  
  for (i in 1:maxiters) {
    beta_last <- beta
    a <- XtWX%*%beta_last
    Dbeta <- as.vector(D%*%beta_last)
    Dbeta_proj <- Proj_sparse(Dbeta, k)
    
    gradient <- as.vector(a - XtWy+rho*t(D)%*%(Dbeta- Dbeta_proj))
    
    Av <- XtWX%*%gradient
    vAv <- t(gradient)%*%Av
    Dv <- as.vector(D%*%gradient)
    
    stepsize <- as.numeric(norm(gradient, "2")^2/(vAv + rho*norm(Dv, "2")^2))
    beta <- beta_last - stepsize*gradient
    
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
    
    rho <- 2*rho
  }
  
  return(list(beta=beta))
}












Proj_sparse <- function(x, k){
  
  res <- sort(abs(x), method="quick", index.return=TRUE)
  ind <- tail(res$ix, k)
  xnew <- rep(0, length(x))
  xnew[ind] <- x[ind]
  return(xnew)
}