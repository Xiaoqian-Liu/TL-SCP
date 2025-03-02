#' Compute a solution path of the all-pairwise penalized estimate for Cox regression. 
#'
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @param X The design matrix.
#' @param penalty The penalty function, options include 'Lasso'(default), 'SCAD', 'MCP', and 'distance'.
#' @param lambdaSeq A decreasing sequence of values for the regularization parameter lambda (for 'Lasso', 'SCAD', and 'MCP').
#' @param nclusters A sequence of number of clusters, only used for penalty = "distance".
#' @param ebic.gm The gamma parameter in EBIC. Default is 0, which corresponds to BIC.
#' @param screening Whether or not to perform feature screening before running the penalized Cox model
#' @param maxiters The maximum number of iterations.
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
allpair_pCox <- function(y, status, X, penalty="Lasso", lambdaSeq=seq(1, 1e-3, length.out=10),
                       nclusters=2:4, ebic.gm = 0, screening = TRUE, maxiters=1e3, tol=1e-4){
  n <- nrow(X)
  p <- ncol(X)
  # if screening is required
  if(screening){
    Y <- survival::Surv(y, status)
    sis_fit <- SIS(X, Y, family='cox', penalty='lasso', tune='bic', varISIS='vanilla')
    nnz_id <- sis_fit$ix
    if(length(nnz_id)>n){
      nnz_num <- floor(p/(4*log(p)))
      sis_fit <- SIS(X, Y, family='cox', penalty='lasso', tune='bic', varISIS='vanilla', nsis = nnz_num)
      nnz_id <- sis_fit$ix
      }
    X0 <- X[, nnz_id] 
  }else{
    X0 <- X
  }
  
  p0 <- ncol(X0)
  D <- allpair_D(p0) 
  
  
  # Note the relation between k in the distance penalty = the number of cluster -1
  fit <- pcox_path(y, status, X0, D, penalty=penalty, lambdaSeq=lambdaSeq,
                   kSeq=nclusters-1, ebic.gm=ebic.gm, maxiters=maxiters, tol=tol)
  
  if(screening){
    Beta_path <- matrix(0, nrow = p, ncol = ncol(fit$Beta))
    for (i in 1:ncol(fit$Beta)) {
      Beta_path[nnz_id, i] <- fit$Beta[, i]
    }
  }else{
    Beta_path <- fit$Beta
  }
  return(list(Beta_path=Beta_path, Iter_path=fit$Iters, IC_path=fit$IC))
}





