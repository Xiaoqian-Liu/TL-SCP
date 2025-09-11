#' Weighted rank averaging
#'
#' @param v1 The first ranking
#' @param v2 The second ranking
#' @param w1 The weight vector associated with v1
#' @param w2 The weight vector associated with v2
#' @export
wRankAveg <- function(v1, v2, w1, w2){
  r1 <- rank(v1, ties.method = "first")
  r2 <- rank(v2, ties.method = "first")

  r <- w1*r1+w2*r2
  tau <- sort(r,index.return=TRUE)$ix
}