#' Compute the solution of weighted LS + fused/generalized Lasso penalty 
#'
#' @param y The vector of response.
#' @param X The design matrix.
#' @param w The vector of weights in the least squares term.
#' @param D The fusion matrix.
#' @param lambda The regularization parameter.
#' @export
fused_Lasso <- function(y, X, w=rep(1, length(y)), D, lambda) {

  n <- length(y)
  p <- ncol(X)
  r <- nrow(D)

  W <- Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X

  
  P <- Matrix(0, nrow = p+2*r, ncol = p+2*r)
  P[1:p, 1:p] <- 0.5*XtWX # should be like this; their code is inconsistent with their paper.
  
  
  q <- c(as.vector(-XtWy), rep(lambda, 2*r))

  A <- cbind(D, Diagonal(n=r, x=-1), Diagonal(n=r, x=1))

  b <- rep(0,r)


  model <- list()

  model$A          <- A
  model$Q          <- P
  model$obj        <- as.vector(q)
  model$modelsense <- 'min'
  model$rhs        <- b
  model$sense      <- rep('=', r)
  model$lb        <- c(rep(-Inf, p), rep(0, 2*r))


  params <- list(OutputFlag=0) 

  result <- gurobi(model, params)

  beta <- result$x[1:p]

  v <- result$x[(p+1):(p+r)]-result$x[(p+r+1):(p+2*r)]
  return(list(beta=beta, v=v))

}
