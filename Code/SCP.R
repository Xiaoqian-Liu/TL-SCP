#' Compute a solution path of the SCP model. 
#'
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @param X The design matrix.
#' @param penalty The penalty function, options include 'Lasso'(default), 'SCAD', 'MCP', and 'distance'.
#' @param lambdaSeq A decreasing sequence of values for the regularization parameter lambda (for 'Lasso', 'SCAD', and 'MCP').
#' @param nclusters A sequence of number of clusters, only used for penalty = "distance".
#' @param ebic.gm The gamma parameter in EBIC. Default is 0, which corresponds to BIC.
#' @param maxiters The maximum number of iterations.
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
SCP <- function(y, status, X, penalty="Lasso", lambdaSeq=seq(1, 1e-3, length.out=10),
                       nclusters=2:4, ebic.gm = 0, maxiters=1e3, tol=1e-4){
  n <- nrow(X)
  p <- ncol(X)
  
  # compute the preliminary estimate using the ridge cox regression
  Y <- cbind(time=y, status=status)
  fit0 <- glmnet(X, Y, family = "cox", alpha = 0, lambda = 1e-4)
  #beta0 <- as.vector(fit0$beta[-1]) # we do not include the gender
  beta0 <- as.vector(fit0$beta) 
  
  # sort X and then use the first-order differencing matrix
  ss <- sort(beta0, index.return = TRUE)
  tau <- ss$ix
  X_order <- X[, tau]
  D <- myGetDkn(k=1, n=ncol(X_order))
  
  
  # Note the relation between k in the distance penalty = the number of cluster -1
  fit <- pcox_path(y, status, X_order, D, penalty=penalty, lambdaSeq=lambdaSeq,
                   kSeq=nclusters-1, ebic.gm=ebic.gm, maxiters=maxiters, tol=tol)
  
  # first reorder the estimates
  beta_path <- fit$Beta
  for (i in 1:ncol(beta_path)) {
    beta_order <- fit$Beta[, i]
    # recover to the original order
    for (j in 1:length(beta_order)) {
      beta_path[tau[j],i] <- beta_order[j]
    }
  }
  
  
  return(list(Beta_path=beta_path, Iter_path=fit$Iters, IC_path=fit$IC))
}





