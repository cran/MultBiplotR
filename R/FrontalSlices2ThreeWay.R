FrontalSlices2ThreeWay <- function (Y, I, J, K){
  # You must provide K, I*J matrices
  Y=as.matrix(Y)
  I=dim(Y)[1]
  JK=dim(Y)[2]
  if (JK != (J*K)) stop("The number of columns of the matrix is not right")
  
  Y3D=array(0, c(I,J,K))
  
  for (k in 1:K){
    Y3D[,,k]=Y[,((k-1)*J+1):(k*J)]
  }
  return(Y3D
)
}
