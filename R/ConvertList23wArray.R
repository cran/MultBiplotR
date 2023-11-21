ConvertList23wArray <- function(X){
  if (!is.list(X)) stop("The input is not a list")
  K=length(X)
  I=dim(X[[1]])[1]
  J=dim(X[[1]])[2]
  
  XX=array(0, c(I,J,K))
  for (k in 1:K){
    if (I!=dim(X[[k]])[1]) stop("Dimensions of all the matrices must be equal")
    if (J!=dim(X[[k]])[2]) stop("Dimensions of all the matrices must be equal")
    XX[,,k]=X[[k]]
  }
  
  dimnames(X)[[3]]=names(XX)
  return(XX)
}