wRankAveg <- function(v1, v2, w1, w2){
  r1 <- rank(v1, ties.method = "first")
  r2 <- rank(v2, ties.method = "first")

  r <- w1*r1+w2*r2
  tau <- sort(r,index.return=TRUE)$ix
}