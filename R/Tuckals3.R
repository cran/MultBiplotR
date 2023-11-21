Tuckals3 <- function(X, Scaling = 5, A=NULL, B=NULL, C=NULL , P=2, Q=2, R=2, tolerance=0.000001, maxiter=1000){
  # Por el momento la transformación inicial será para las variables que, generalmente, se colocan en el segundo modo
  # Transformaremos cada variable en cada ocasión, es decir, para cada combinación variable-ocasión centramos sobre los individuos
  # Dependiendo de la transformación los biplots obtenidos pueden ser distintos
  
  dims=dim(X)
  I=dims[1]
  J=dims[2]
  K=dims[3]
  Inames=dimnames(X)[[1]]
  Jnames=dimnames(X)[[2]]
  Knames=dimnames(X)[[3]]

  X1=ThreeWay2FrontalSlices(X, 1)$FrontalSlice
  X2=ThreeWay2FrontalSlices(X, 2)$FrontalSlice
  X3=ThreeWay2FrontalSlices(X, 3)$FrontalSlice
  
  # Initial values for A, B, C y G
  if (is.null(A)) A=svd(X1)$u[,1:P]
  if (is.null(B)) B=svd(X2)$u[,1:Q]
  if (is.null(C)) C=svd(X3)$u[,1:R]
  
  K1=kronecker(C,B)
  G=t(A) %*% X1 %*% K1
  Xest=A %*% G %*% t(K1)
  
  SCRold = sum((X1-Xest)^2)
  error=1
  iter=0
  while ((error>tolerance) & (iter<maxiter)){
    iter=iter+1
    A=svd(X1 %*% K1)$u[,1:P]
    K2=kronecker(C,A)
    B=svd(X2 %*% K2)$u[,1:Q]
    K3=kronecker(A, B)
    C=svd(X3 %*% K3)$u[,1:R]
    K1=kronecker(C,B)
    G=t(A) %*% X1 %*% K1
    Xest=A %*% G %*% t(K1)
    SCRnew = sum((X1-Xest)^2)
    error=abs(SCRold-SCRnew)
    SCRold=SCRnew
    cat("Iteration :", iter, " -  Error:",error, "\n")
  }
  
  SCE=sum(G^2)
  SCT=sum(X1^2)
  ExPer=100*SCE/SCT
  
  rownames(A)=Inames
  colnames(A)=paste("Dim:",1:P, sep="")
  
  rownames(B)=Jnames
  colnames(B)=paste("Dim:",1:Q, sep="")
  
  rownames(C)=Knames
  colnames(C)=paste("Dim:",1:R, sep="")
  # Devolvemos también el biplot interactivo dejando el modo 1 en filas
  
  BiplotY = list()
  BiplotY$Title = "Interactive Biplot Tuckals 3"
  BiplotY$Type = "Interactive Biplot" 
  BiplotY$Non_Scaled_Data = X1
  BiplotY$alpha=0
  BiplotY$Dimension=P
  BiplotY$Means = apply(X1, 2, mean)
  BiplotY$Medians = apply(X1, 2, median)
  BiplotY$Deviations = apply(X1, 2, sd)
  BiplotY$Minima = apply(X1, 2, min)
  BiplotY$Maxima = apply(X1, 2, max)
  BiplotY$P25 = apply(X1, 2, quantile)[2, ]
  BiplotY$P75 = apply(X1, 2, quantile)[4, ]
  BiplotY$GMean = mean(X1)
  BiplotY$Initial_Transformation = Scaling
  Data = InitialTransform(X1, transform = Scaling)
  Y = Data$X
  BiplotY$Scaled_Data = Y
  
  BiplotY$nrows = I
  BiplotY$ncols = J*K
  BiplotY$dim = P
  
  Y=BiplotY$Scaled_Data
  sct=sum(Y^2)
  
  scf = apply((Y^2), 1, sum)
  scc = apply((Y^2), 2, sum)
  
  a=A
  b=t(G %*% t(K1))
  sce=rep(0,P)
  cf=matrix(0,I,P)
  cc=matrix(0,J*K,P)
  for (i in 1:P){
    Yesp=matrix(a[,i])%*%t(matrix(b[,i]))
    scfe = apply((Yesp^2), 1, sum)
    scce = apply((Yesp^2), 2, sum)
    scfr = apply((Y-Yesp^2), 1, sum)
    sccr = apply((Y-Yesp^2), 2, sum)
    cf[,i]=scfe/scf
    cc[,i]=scce/scc
    sce[i]=sum(Yesp^2)
  }
  
  cf=100*cf
  cc=100*cc
  BiplotY$EigenValues = sce
  BiplotY$Inertia = 100 * sce/sct
  
  BiplotY$CumInertia = cumsum(BiplotY$Inertia)
  
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/I
  scb = scb/J*K
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf
  
  rownames(b) = colnames(X1)
  colnames(a) = paste("Dim", 1:P)
  
  BiplotY$RowCoordinates=a
  BiplotY$ColCoordinates=b
  
  rownames(cf) = rownames(X1)
  colnames(cf) = paste("Dim", 1:P)
  BiplotY$RowContributions=cf
  
  rownames(cc) = colnames(X1)
  colnames(cc) = paste("Dim", 1:P)
  BiplotY$ColContributions=cc
  
  BiplotY$Scale_Factor = scf
  BiplotY$ClusterType="us"
  BiplotY$Clusters = as.factor(matrix(1,nrow(BiplotY$RowContributions), 1))
  BiplotY$ClusterColors="blue"
  BiplotY$ClusterNames="Cluster 1"
  
  class(BiplotY)="ContinuousBiplot"
  
  Result=list(X=X, Xest=Xest, A=A, B=B, C=C, G=G, RSS=SCRnew, ESS=SCE, TSS=SCT, PercentExplained=ExPer,
              Biplot=BiplotY)
  class(Result)="Tuckals3"
  return(Result)
}



