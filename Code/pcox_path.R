#' Fit a solution path for the fusion penalized Cox regression model.
#'
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @param X The design matrix.
#' @param D The fusion matrix in the penalty term.
#' @param penalty The penalty function, options include 'Lasso'(default), 'SCAD', 'MCP', and 'distance'.
#' @param lambdaSeq A decreasing sequence of values for the regularization parameter lambda (for 'Lasso', 'SCAD', and 'MCP').
#' @param kSeq A sequence of increasing values of k for  penalty = "distance".
#' @param ebic.gm The gamma parameter in EBIC. Default is 0, which corresponds to BIC. 
#' @param maxiters The maximum number of iterations.
#' @param tol The tolerance value for early stopping the algorithm. 
#' @export
pcox_path <- function(y, status, X, D, penalty="Lasso", lambdaSeq=seq(1, 1e-3, length.out=10),
                      kSeq=1:3, ebic.gm = 0, maxiters=1e3, tol=1e-4){
  
  n <- nrow(X)
  p <- ncol(X)
  
  #======different penalty options
  ####### Lasso or SCAD or MCP
  if(penalty == 'Lasso' || penalty == 'SCAD' || penalty == 'MCP'){
    nlam <- length(lambdaSeq)
    # setup the output
    Beta <- matrix(NA, nrow = p, ncol = nlam)
    Iters <- rep(NA, nlam)
    IC <- rep(NA, nlam)

    Y <- cbind(time=y, status=status)
    fit_cox <- glmnet(X, Y, family = "cox", alpha = 0, lambda = 1e-4)
    beta0 <- runmed(fit_cox$beta, k=5) #median filter
    #beta0 <- fit_cox$beta
    
    for(i in 1:nlam){
      lambda <- lambdaSeq[i]
      fit <- pcox(y, status, X, D, penalty=penalty, lambda = lambda, beta0 = beta0, maxiters = maxiters, tol = tol)
      # save the solution and update the initial value
      beta0 <- Beta[, i] <- fit$beta
      # save number of iterations
      Iters[i] <- fit$iters
      # compute the EBIC
      nu <- length(which(abs(D%*%beta0)>1e-3))
      tau <- choose(nrow(D), nu)
      IC[i] <- -2*pll_cox(y, status, X, beta0) + log(n)*nu + 2*ebic.gm*log(tau)
    }
  }
  
  
  
  ####### Distance
  if(penalty == 'distance'){
    nk <- length(kSeq)
    # setup the output
    Beta <- matrix(NA, nrow = p, ncol = nk)
    Iters <- rep(NA, nk)
    IC <- rep(NA, nk)
    # Needs good initial values for distance penalization
    Y <- cbind(time=y, status=status)
    fit_cox <- glmnet(X, Y, family = "cox", alpha = 0, lambda = 1e-4)
    beta0 <- runmed(fit_cox$beta, k=5) #median filter
    beta0 <- as.numeric(as.matrix(fit_cox$beta))
    #beta0 <- fit_cox$beta
    
    for(i in 1:nk){
      fit <- pcox(y, status, X, D, penalty="distance", k=kSeq[i], beta0 = beta0, maxiters = maxiters, tol = tol)
      # save the solution
      Beta[, i] <- as.vector(fit$beta)
      # save number of iterations
      Iters[i] <- fit$iters
      # compute the EBIC
      nu <- length(which(abs(D%*%fit$beta)>1e-3))
      tau <- choose(nrow(D), nu)
      IC[i] <- -2*pll_cox(y, status, X, fit$beta) + log(n)*nu + 2*ebic.gm* log(tau)
    }
  }
  
  # return the results
  return(list(Beta=Beta, Iters=Iters, IC=IC))
}
