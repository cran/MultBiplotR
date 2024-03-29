binBvn <- function(rho, row.cuts, col.cuts, bins=4){
  row.cuts <- if (missing(row.cuts)) c(-Inf, 1:(bins - 1)/bins, Inf) else  c(-Inf, row.cuts, Inf)
  col.cuts <- if (missing(col.cuts)) c(-Inf, 1:(bins - 1)/bins, Inf) else  c(-Inf, col.cuts, Inf)
  r <- length(row.cuts) - 1
  c <- length(col.cuts) - 1
  P <- matrix(0, r, c)
  R <- matrix(c(1, rho, rho, 1), 2, 2)
  for (i in 1:r){
    for (j in 1:c){
      P[i,j] <- pmvnorm(lower=c(row.cuts[i], col.cuts[j]),
                                 upper=c(row.cuts[i+1], col.cuts[j+1]),
                                 corr=R)
    }
  }
  P
}
