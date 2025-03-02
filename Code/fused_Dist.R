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
fused_iter_Dist <- function(y, X, w=rep(1, length(y)), D, k, rho=1,  beta0=rep(0, ncol(X)), maxiters=1e2, tol=1e-5){
  n <- length(w)
  W <- Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X
  beta <- beta0
  
  for (i in 1:maxiters) {
    # update some terms
    XtW <- t(X)%*%W
    XtWy <- XtW%*%y
    XtWX <- XtW%*%X
    
    
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
    
    # update X and D based on the new beta
    ss <- sort(beta, index.return = TRUE)
    tau <- ss$ix
    X <- X[, tau]
    D <- myGetDkn(k=1, n=ncol(X))
    
  }
  
  return(list(beta=beta))
}







#' Compute the solution of weighted LS + double fusion-distance penalties 
#'
#' @param y The vector of response
#' @param X The design matrix
#' @param w The weights in the least squares term
#' @param D_in One fusion matrix (for the within-segment penalty in aCARDS)
#' @param D_bt The other fusion matrix (for the between-segment penalty in aCARDS)
#' @param rho The regularization parameter, typically a large value
#' @param k_in The constraint parameter (number of nonzeros) in the within-segment penalty. 
#' @param k_bt The constraint parameter (number of nonzeros) in the between-segment penalty. 
#' @param beta0 The initial vector of coefficients
#' @param maxiters Maximum number of iterations
#' @param tol The tolerance value to terminate the algorithm 
#' @export
fused_2Dist <- function(y, X, w, D_in, D_bt, rho, k_in, k_bt, beta0, maxiters=1e2, tol=1e-6){
  n <- length(w)
  W <- Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X
  beta <- beta0
  
  for (i in 1:maxiters) {
    beta_last <- beta
    a <- XtWX%*%beta_last
    Dbeta_in <- as.vector(D_in%*%beta_last)
    Dbeta_bt <- as.vector(D_bt%*%beta_last)
    Dbeta_in_proj <- Proj_sparse(Dbeta_in, k_in)
    Dbeta_bt_proj <- Proj_sparse(Dbeta_bt, k_bt)
    
    Dbeta <- c(Dbeta_in, Dbeta_bt)
    Dbeta_proj <- c(Dbeta_in_proj, Dbeta_bt_proj)
    D <- rbind(D_in, D_bt)
    
    gradient <- as.vector(a - XtWy+rho*t(D)%*%(Dbeta- Dbeta_proj))
    
    Av <- XtWX%*%gradient
    vAv <- t(gradient)%*%Av
    Dv <- as.vector(D%*%gradient)
    
    stepsize <- as.numeric(norm(gradient, "2")^2/(vAv + rho*norm(Dv, "2")^2))
    beta <- beta_last - stepsize*gradient
    
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
    
    rho <- 2*rho
  }
  
  return(list(beta=beta, iters=i))
}










Proj_sparse <- function(x, k){
  
  res <- sort(abs(x), method="quick", index.return=TRUE)
  ind <- tail(res$ix, k)
  xnew <- rep(0, length(x))
  xnew[ind] <- x[ind]
  return(xnew)
}