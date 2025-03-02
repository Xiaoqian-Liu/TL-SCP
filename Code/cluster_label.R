#' Create cluster labels for a given vector of coefficients in a data-driven manner.
#'
#' @param beta The vector of coefficients.
#' @export
cluster_label <- function(beta){
  threshold_value <- max(abs(beta))/100
  ss <- sort(beta, index.return = TRUE)
  sorted_beta <- ss$x
  #threshold_value<- quantile(diff(sorted_beta), 0.9)
  ind_set <- which(diff(sorted_beta) > threshold_value)
  label_vec <- rep(NA, length(beta))
  ind_set <- c(0, ind_set, length(beta))
  for (i in 1:(length(ind_set)-1)) {
    sorted_idx_set <- (ind_set[i]+1):(ind_set[i+1])
    original_idx_set <- ss$ix[sorted_idx_set]
    label_vec[original_idx_set] <- i
  }
  return(label_vec)
}




#' Create cluster labels for a given vector of coefficients from the distance penalization
#'
#' @param beta The vector of coefficients.
#' @param k The number of clusters
#' @export
cluster_label_dist <- function(beta, k){
  
  ss <- sort(beta, index.return = TRUE)
  sorted_beta <- ss$x
  top_diff <- head(sort(diff(sorted_beta), decreasing=TRUE), k-1)
  top_loc <- sort(match(top_diff, diff(sorted_beta)))
 
  label_vec <- rep(NA, length(beta))
  ind_set <- c(0, top_loc, length(beta))
  for (i in 1:(length(ind_set)-1)) {
    sorted_idx_set <- (ind_set[i]+1):(ind_set[i+1])
    original_idx_set <- ss$ix[sorted_idx_set]
    label_vec[original_idx_set] <- i
  }
  return(label_vec)
}





#' Create cluster labels for a given vector of coefficients in a data-driven manner.
#'
#' @param beta The vector of coefficients.
#' @param k The number of clusters
#' @export
cluster_label_kmeans <- function(beta, k){
  
  if(length(k)==1){
    kmeans_results <- kmeans(beta, k) 
    label_vec <- kmeans_results$cluster
  }else{
    numK <- length(k)
    loss <- Inf
    for (j in 1:numK) {
      kmeans_results <- kmeans(beta, k[j]) 
      if(kmeans_results$totss<loss){
        label_vec <- kmeans_results$cluster
      }
    }
  }

  
  return(label_vec)
}

