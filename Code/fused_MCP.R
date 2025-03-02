#' Compute the solution of weighted LS + fused MCP penalty. 
#'
#' @param y The response vector.
#' @param X The design matrix.
#' @param w The vector of weights in the least squares term.
#' @param D The fusion matrix.
#' @param lambda The regularization parameter.
#' @param a The hyper-parameter in MCP.
#' @param beta0 The initial vector of coefficients.
#' @param maxiters The maximum number of iterations for the algorithm
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
fused_MCP <- function(y, X, w, D, lambda, a = 3, beta0 = rep(0, ncol(X)), maxiters=1e3, tol=1e-5) {
  
  # create a vector to record current stage beta
  cur_beta <- beta0
  
  # then initialize the algorithm
  for (i in 1:maxiters) {
    
    # obtain fused beta
    Dbeta <- D %*% cur_beta
    
    # obtain the diagonal prime matrix
    prime_mat <- diag(as.vector(MCP_prime(Dbeta, lambda, a)))
    
    # construct the new D matrix 
    new_D <- prime_mat %*% D
    
    # obtain next stage beta by fused lasso with lambda = 1
    next_beta <- fused_Lasso(y, X, w, D = new_D, lambda = 1)$beta
  
    # check convergence 
    A <- norm(next_beta-cur_beta, "2")/(1+norm(cur_beta, "2"))
    
    if(A<tol){
      # terminate the loop
      break
    }
    
    # update for the next iteration
    cur_beta <- next_beta
    
  }
  
  return(list(beta=cur_beta, iters=i))
  
}



#' Compute the solution of weighted LS + double fused MCP penalties.
#'
#' @param y The response vector.
#' @param X The design matrix.
#' @param w The vector of weights in the least squares term.
#' @param D_in One fusion matrix (for the within-segment penalty in aCARDS).
#' @param D_bt The other fusion matrix (for the between-segment penalty in aCARDS).
#' @param lambda_in The regularization parameter associated with D_in.
#' @param lambda_bt The regularization parameter associated with D_bt.
#' @param a The hyper-parameter in MCP.
#' @param beta0 The initial vector of coefficients.
#' @param maxiters The maximum number of iterations for the algorithm.
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
fused_2MCP <- function(y, X, w=rep(1, length(y)), D_in, D_bt, lambda_in, lambda_bt, 
                        a=3.7, beta0=rep(0, ncol(X)), maxiters=1e2, tol=1e-5) {
  
  # the key is that we transform it to a fused lasso problem
  
  for (k in 1:maxiters) {
    
    # for each iteration, D changes
    Dbeta_in <- D_in%*%beta0
    Dbeta_bt <- D_bt%*%beta0
    
    
    # obtain the two L diagonal matrices, whose diagonal values are MCP derivative
    L_in <- Diagonal(n=nrow(D_in), x=MCP_prime(Dbeta_in, lambda=lambda_in, a=a))
    L_bt <- Diagonal(n=nrow(D_bt), x=MCP_prime(Dbeta_bt, lambda=lambda_bt, a=a))
    D1 <- L_in%*%D_in
    D2 <- L_bt%*%D_bt
    D <- rbind(D1, D2)
    
    # at each iteration, solve a fused lasso problem
    beta <- fused_Lasso(y, X, w, D=D, lambda=1)$beta #Note here lambda=1 for the fused lasso problem
    
    #check convergence
    if ( norm(beta-beta0,'2') < tol*(1 + norm(beta0,'2')) ) break
    
    # update the estimate
    beta0 <- beta
  }
  
  return(list(beta=beta, iters=k))
}











#' Compute the derivative of MCP. 
#'
#' @param beta The vector of coefficients.
#' @param lambda The regularization parameter.
#' @param a The hyper-parameter in MCP.
#' @export
MCP_prime <- function(beta, lambda, a){
  
  # get the beta size
  size = length(beta)
  
  # create a row vector to store the first derivative for each beta entries
  prime_result = matrix(0, nrow = 1, ncol = size)
  
  # the MCP derivative threshold
  thresh = lambda * a
  
  for (i in 1:size) {
    
    prime_result[1, i] = ifelse(abs(beta[i, 1]) <= thresh, (lambda - abs(beta[i, 1])/a) * sign(beta[i, 1]), 0)  
    
  }
  
  return(prime_result)
  
}

