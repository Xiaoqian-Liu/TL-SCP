#' Compute the solution path for transfer learning SCP. 
#'
#' @param Y_T The survial info of the target dataset, including the survival time in the first column
#'            and the status info in the second column. 
#' @param X_T The design matrix of the target dataset. 
#' @param Y_S The survial info of the source dataset, including the survival time in the first column
#'            and the status info in the second column. 
#' @param X_S The design matrix of the source dataset. 
#' @param penalty The penalty function, options include 'Lasso'(default), 'SCAD', 'MCP', and 'distance'.
#' @param lambdaSeq A decreasing sequence of values for the regularization parameter lambda (for 'Lasso', 'SCAD', and 'MCP').
#' @param nclusters A sequence of number of clusters, only used for penalty = "distance".
#' @param ebic.gm The gamma parameter in EBIC. Default is 0, which corresponds to BIC.
#' @param maxiters The maximum number of iterations.
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
trans_SCP <- function(Y_T, X_T, Y_S, X_S, penalty="Lasso", lambdaSeq=seq(1, 1e-3, length.out=10),
                       nclusters=2:4, ebic.gm = 0, maxiters=1e3, tol=1e-4){

  n_S <- nrow(X_S)
  n_T <- nrow(X_T)
  n <- n_S+n_T
  
  p <- ncol(X_T)
  
  # compute the preliminary estimate using the source data
  fit_S <- glmnet(X_S, Y_S, family = "cox", alpha = 0, lambda = 0)
  beta_S <- as.vector(fit_S$beta[1:p]) 
  
  # compute the preliminary estimate using the target data
  fit_T <- glmnet(X_T, Y_T, family = "cox", alpha = 0, lambda = 0)
  beta_T <- as.vector(fit_T$beta) 
  
  tau <- wRankAveg(beta_S, beta_T, n_S/n, n_T/n)
  
  X_order <- X_T[, tau]
  D <- myGetDkn(k=1, n=ncol(X_order))
  
  
  # Note the relation between k in the distance penalty = the number of cluster -1
  fit <- pcox_path(y=Y_T[, 1], status=Y_T[, 2], X_order, D, penalty=penalty, lambdaSeq=lambdaSeq,
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





