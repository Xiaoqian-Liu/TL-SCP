#' Compute all-pairwise differening matrix.
#'
#' @param n The number of predictors.
#' @export
allpair_D <- function(n) {
  A <- combn(n, 2) # get all pairs
  m <- ncol(A)
  D <- Matrix(0, nrow = m, ncol = n, sparse = TRUE)
  for (i in 1:m) {
    D[i, A[, i]] <- c(-1, 1)
  }
  return(D)
}
