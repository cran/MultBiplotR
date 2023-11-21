Convert3wArray2List= function(X){
  if (length(dim(X))!=3) stop("The input must be a 3-way array")
  I=dim(X)[1]
  J=dim(X)[2]
  K=dim(X)[3]
  
  XX=list()
  for (k in 1:K){
    XX[[k]]=X[,,k]
  }
  names(XX)=dimnames(X)[[3]]
  return(XX)

}

