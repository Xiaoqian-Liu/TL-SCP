#' Compute w and z in the transformed weighted LS problem for Cox regression
#'
#' @param eta The eta vector in the definition of w and z, eta=Xbeta.
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @export
update_wz <- function(eta, y, status){
  
  # get the number of sample size, which is the dimension of w and z
  n <- length(y)
  w <- rep(0, n)
  z <- rep(NA, n)
  
  for (k in 1:n) {
    # compute C_k, the set of i with t_i < y_k
    indset <- which(y<=y[k] & status==1)
    
    if(length(indset)==0){
      
      w[k] <- 0 # before is 1.
      z[k] <- eta[k] #+ (status[k] - 1)
      
    }else{
      # compute w and a part in z
      s <- 0
      for(i in indset){
        # compute R_i, the set of indices j with y_j >= t_i
        R <- which(y>= y[i])
        c <- sum(exp(eta[R]))
        a <- exp(eta[k])
        if(is.infinite(c) || is.na((a*c - a^2)/c^2)){
          w[k] <- w[k]+0
        }else{
          w[k] <- w[k]+ (a*c - a^2)/c^2
        }

        
        s <- s+a/c
      }
      
      # compute z
      z[k] <- eta[k] + (status[k] - s)/w[k]
    }
    
    
  }
  
  return(list(w=w, z=z))
  
}