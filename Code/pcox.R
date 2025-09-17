#' Fit the penalized Cox regression model at a fixed value of the penalization parameter. 
#'
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @param X The design matrix.
#' @param D The fusion matrix in the penalty term. 
#' @param penalty The penalty function, options include 'Lasso'(default), 'SCAD', 'MCP', and 'distance'.
#' @param lambda The regularization parameter lambda.
#' @param k The regularization parameter k for the distance penalty (equal to the number of nonzeros in Dbeta).
#' @param beta0 The initial vector of coefficients. 
#' @param maxiters The maximum number of iterations.
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
pcox <- function(y, status, X, D, penalty="Lasso", lambda=0, k=2, rho=1e-3, beta0=rep(0, ncol(X)), maxiters=1e3, tol=1e-4){
  
  n <- length(y)
  # get the initial value of eta
  eta0 <- X%*%beta0
  
  for(i in 1:maxiters){
    # compute the Hessian matrix W and z through a function
    wz <- update_wz(eta0, y, status)
    w <- wz$w
    z <- wz$z
    
    
    # update beta through a fused penalized LS problem
    if(penalty == 'Lasso'){
      beta <- fused_Lasso(y=z, X=X, w=w/n, D=D, lambda=lambda)$beta
    }
    
    if(penalty == 'SCAD'){
      beta <- fused_SCAD(y=z, X=X, w=w/n, D=D, lambda=lambda, beta0=beta0, maxiters=maxiters, tol=tol)$beta
    }
    
    if(penalty == 'MCP'){
      beta <- fused_MCP(y=z, X=X, w=w/n, D=D, lambda=lambda, beta0=beta0, maxiters=maxiters, tol=tol)$beta
    }
    
    if(penalty == 'distance'){
      beta <- fused_Dist(y=z, X=X, w=w/n, D=D, k=k, rho=rho, beta0=beta0, maxiters=maxiters, tol=tol)$beta
    }
   
    
    
    # check convergence 
    A <- norm(beta-beta0, "2")/(1+norm(beta0, "2"))
    if(A<tol){
      # terminate the loop
      break
    }else{
      # update for the next iteration
      beta0 <- beta
      eta0 <- X%*%beta0
    }
  }
 
  # return the results
  return(list(beta=beta, iters=i))
}





