ThreeWay2FrontalSlices <- function(X, Slice=1){
  # Index es el modo que debe colocarse en Fila de la matriz.
  # Los otros dos se codifican Interactivamente
  
  dims=dim(X)
  I=dims[1]
  J=dims[2]
  K=dims[3]
  
  Inames=dimnames(X)[[1]]
  Jnames=dimnames(X)[[2]]
  Knames=dimnames(X)[[3]]
  
  if (Slice==1){
    FS = X[,1,]
    colnames(FS)=paste(Jnames[1], "-", Knames,  sep="")
    for (j in 2:J){
      XX=X[,j,] 
      colnames(XX)=paste(Jnames[j], "-", Knames,  sep="")
      FS = cbind(FS, XX)
    }
  }
  
  if (Slice==2){
    FS = X[1,,]
    colnames(FS)=paste(Inames[1], "-", Knames,  sep="")
    for (i in 2:I){
      XX=X[i,,]
      colnames(XX)=paste(Inames[i], "-", Knames,  sep="")
      FS = cbind(FS, XX)
    }
  }
  
  if (Slice==3){
    FS = t(X[1,,])
    colnames(FS)=paste(Inames[1], "-", Jnames,  sep="")
    for (i in 2:I){
      XX=t(X[i,,])
      colnames(XX)=paste(Inames[i], "-", Jnames,  sep="")
      FS = cbind(FS, XX)
    }
  }
  
  return(list(FrontalSlice=FS, I=I, J=J, K=K, Slice=Slice))
  
}
