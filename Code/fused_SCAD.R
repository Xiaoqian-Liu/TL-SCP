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




#' Compute the solution of weighted LS + double fused SCAD penalties.
#'
#' @param y The response vector.
#' @param X The design matrix.
#' @param w The vector of weights in the least squares term.
#' @param D_in One fusion matrix (for the within-segment penalty in aCARDS).
#' @param D_bt The other fusion matrix (for the between-segment penalty in aCARDS).
#' @param lambda_in The regularization parameter associated with D_in.
#' @param lambda_bt The regularization parameter associated with D_bt.
#' @param a The hyper-parameter in SCAD.
#' @param beta0 The initial vector of coefficients.
#' @param maxiters The maximum number of iterations for the algorithm.
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
fused_2SCAD <- function(y, X, w=rep(1, length(y)), D_in, D_bt, lambda_in, lambda_bt, 
                        a=3.7, beta0=rep(0, ncol(X)), maxiters=1e2, tol=1e-5) {
  # the key is that we transform it to a fused lasso problem
  
  for (k in 1:maxiters) {
    # for each iteration, D changes
    Dbeta_in <- D_in%*%beta0
    Dbeta_bt <- D_bt%*%beta0
    L_in <- Diagonal(n=nrow(D_in), x=scad_prime(Dbeta_in, lambda=lambda_in, a=a))
    L_bt <- Diagonal(n=nrow(D_bt), x=scad_prime(Dbeta_bt, lambda=lambda_bt, a=a))
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




# set.seed(123)
# n <- 100
# p <- 10
# beta <- rep(0, p)
# beta[1:3] <- 1:3
# 
# X <- matrix(data = rnorm(n*p), nrow = n)
# y <- X%*%beta + rnorm(n)
# 
# # test Lasso, gurobi is better than OSQP
# D <- Diagonal(n=p, x=1)
# fit_lasso<- sol_fused_Lasso(y, X, w=rep(1, ncol(X)), D, lambda=0.1)
# round(fit_lasso$beta, 8)
# fit1_lasso <- glmnet(X, y, weights=rep(1/sqrt(n), length(y)), alpha = 1, lambda = 0.1, intercept = FALSE)
# as.vector(fit1_lasso$beta)
# 
# 
# 
# library(ncvreg)
# fit1 <- ncvreg(X, y, penalty = 'SCAD', lambda = 0.1)
# as.vector(fit1$beta[-1])
# 
# D <- Diagonal(n=p, x=1)
# fit2 <- fused_SCAD(y, X, w=rep(1, ncol(X)), D=D, lambda = 0.1, tol = 1e-8)
# 
# round(fit2$beta, 8)
# as.vector(fit1$beta[-1])
