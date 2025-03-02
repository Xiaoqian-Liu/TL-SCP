#' Compute the partial log-likelihood for Cox regression.
#'
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @param X The design matrix.
#' @param beta The vector of coefficients. 
#' @param weight The vector of weights associate with each observation, default is 1.
#' @export
pll_cox <- function(y, status, X, beta, weight = rep(1, length(y))){
  fail_times <- unique(y[status==1])
  m <- length(fail_times)
  # Define a vector to record likelihood at each specific failure time
  likelihood_vec <- rep(NA, m)

  for(j in 1:m){
    D <- which(y==fail_times[j] & status==1) # consider ties
    R <- which(y>=fail_times[j])
    # Compute the current stage likelihood with the given weights
    Xb1 <- X[D, ] %*% beta
    part1 <- sum(diag(weight[D], nrow = length(D)) %*% Xb1)
    Xb2 <- X[R, ] %*% beta
    part2 <- sum(weight[D]) * log_sum_exp_weighted(x = Xb2, w = weight[R]) 

    likelihood_vec[j] = part1 - part2
  }

  # Return the partial log-likelihood
  return(sum(likelihood_vec))
}


#' A function to apply the log-sum-exp trick
#'
#' @param X the parameter we wants to compute the log exp sum
#' @param w The weights of the exp
#' @export
log_sum_exp_weighted <- function(x, w) {
  a <- max(x)  # Find the maximum value to scale the computation
  sum_exp <- sum(w * exp(x - a))  # Compute the weighted sum of exponentials in a numerically stable way
  result <- a + log(sum_exp)  # Adjust back by adding 'a' and applying log
  return(result)
}






# old function

# pll_cox <- function(y, status, X, beta, weight = rep(1, length(y))){
#   fail_times <- unique(y[status==1])
#   m <- length(fail_times)
#   # Define a vector to record likelihood at each specific failure time
#   likelihood_vec <- rep(NA, m)
# 
#   for(j in 1:m){
#     D <- which(y==fail_times[j] & status==1) # consider ties
#     R <- which(y>=fail_times[j])
#     # Compute the current stage likelihood with the given weights
#     Xb1 <- X[D, ] %*% beta
#     numerator <- exp(sum(diag(weight[D], nrow = length(D)) %*% Xb1))
#     Xb2 <- X[R, ] %*% beta
#     denominator <- sum(diag(weight[R], nrow = length(R)) %*% exp(Xb2))^sum(weight[D])
#     
#     print(diag(weight[R], nrow = length(R)) %*% exp(Xb2))
#     likelihood_vec[j] = numerator/denominator
#   }
# 
# 
# 
#   # Return the partial log-likelihood
#   return(log(prod(likelihood_vec)))
# }


