Tucker3PCovR <- function(X , Y3D, J=NULL, K=NULL, r1=2 , r2=2 , r3=2 , Input="Three",
                         conv=1e-6 , OriginalAlfa = 0.5, AlternativeLossF=1 ,
                         nRuns=100 , StartSeed=0, ScalingX = 5, ScalingY = 5){
  
  if (Input=="Three"){
    J=dim(Y3D)[2]
    K=dim(Y3D)[3]
  }
  else{
    if (is.null(J) | is.null(K)) stop("You must provide the dimensions of the ThreeWay array")
  }
  
  I=dim(X)[1]
  L=dim(X)[2]
  
  if (is.null(rownames(X))) {rownames(X)=paste("i",1:I,sep="")}
  if (is.null(colnames(X))) {colnames(X)=paste("X",1:L,sep="")}
  
  
  
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering",
                              "Column centering", "Standardize columns", "Row centering",
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
  
  if (is.numeric(ScalingX))
    ScalingX = ContinuousDataTransform[ScalingX]
  if (is.numeric(ScalingY))
    ScalingY = ContinuousDataTransform[ScalingY]
  
  X = as.matrix(X)
  
  # The initial matrices must be transformed according to the Scaling parameters.
  BiplotX = list()
  BiplotX$Title = "T3PcovR Biplot"
  BiplotX$Type = "T3PcovR Biplot" 
  BiplotX$Non_Scaled_Data = X
  BiplotX$alpha=0
  BiplotX$Dimension=r1
  BiplotX$Means = apply(X, 2, mean)
  BiplotX$Medians = apply(X, 2, median)
  BiplotX$Deviations = apply(X, 2, sd)
  BiplotX$Minima = apply(X, 2, min)
  BiplotX$Maxima = apply(X, 2, max)
  BiplotX$P25 = apply(X, 2, quantile)[2, ]
  BiplotX$P75 = apply(X, 2, quantile)[4, ]
  BiplotX$GMean = mean(X)
  BiplotX$Initial_Transformation = ScalingX
  Data = InitialTransform(X, transform = ScalingX)
  X = Data$X
  BiplotX$Scaled_Data = X
  class(BiplotX)="ContinuousBiplot"
  rm(Data)
  

  Y = as.matrix(Y3D)
  if (Input=="Three") Y = matrix(Y , I , J * K) # P slices are concatenated horizontally
  
  if (is.null(rownames(Y))) {rownames(Y)=rownames(X)}
  if (is.null(colnames(Y))) {
    names=NULL
    for (k in 1:K)
      names=c(names, paste("V",1:J,"-O",k, sep=""))
    colnames(Y)=names
      }
  
  BiplotY = list()
  BiplotY$Title = "T3PcovR Biplot"
  BiplotY$Type = "T3PcovR Biplot" 
  BiplotY$Non_Scaled_Data = Y
  BiplotY$alpha=0
  BiplotY$Dimension=min(r1, r2, r3)
  BiplotY$Means = apply(Y, 2, mean)
  BiplotY$Medians = apply(Y, 2, median)
  BiplotY$Deviations = apply(Y, 2, sd)
  BiplotY$Minima = apply(Y, 2, min)
  BiplotY$Maxima = apply(Y, 2, max)
  BiplotY$P25 = apply(Y, 2, quantile)[2, ]
  BiplotY$P75 = apply(Y, 2, quantile)[4, ]
  BiplotY$GMean = mean(Y)
  BiplotY$Initial_Transformation = ScalingY
  Data = InitialTransform(Y, transform = ScalingY)
  Y = Data$X
  BiplotY$Scaled_Data = Y

  
  t3pc=t3pcovr(X=X , Y=Y , I=I , J=J , K =K, L=L , r1=r1 , r2=r2 , r3=r3 ,
                       conv = conv, OriginalAlfa = 0.5, AlternativeLossF=1 ,
                       nRuns=100 , StartSeed=0)

  # Construcción del biplot para X (separado)
  
  BiplotX$nrows = I
  BiplotX$ncols = L
  BiplotX$dim = r1
  
  X=BiplotX$Scaled_Data
  sct=sum(X^2)
  
  scf = apply((X^2), 1, sum)
  scc = apply((X^2), 2, sum)
  
  a=t3pc$A
  rownames(a)=rownames(X)
  b=t3pc$B2
  sce=rep(0,r1)
  cf=matrix(0,I,r1)
  cc=matrix(0,L,r1)
  for (i in 1:r1){
    Xesp=matrix(a[,i])%*%t(matrix(b[,i]))
    scfe = apply((Xesp^2), 1, sum)
    scce = apply((Xesp^2), 2, sum)
    cf[,i]=scfe/scf
    cc[,i]=scce/scc
    sce[i]=sum(Xesp^2)
  }
  cf=100*cf
  cc=100*cc
  BiplotX$EigenValues = sce
  BiplotX$Inertia = 100 * sce/sct
  
  BiplotX$CumInertia = cumsum(BiplotX$Inertia)
  BiplotX$EV = t3pc$B2
  
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/I
  scb = scb/L
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf
  
  BiplotX$RowCoordinates=a
  BiplotX$ColCoordinates=b
  

  rownames(cf) = rownames(X)
  colnames(cf) = paste("Dim", 1:r1)
  BiplotX$RowContributions=cf
  
  rownames(cc) = colnames(X)
  colnames(cc) = paste("Dim", 1:r1)
  BiplotX$ColContributions=cc
  
  BiplotX$Scale_Factor = scf
  BiplotX$ClusterType="us"
  BiplotX$Clusters = as.factor(matrix(1,nrow(BiplotX$RowContributions), 1))
  BiplotX$ClusterColors="blue"
  BiplotX$ClusterNames="Cluster 1"
  
  class(BiplotX)="ContinuousBiplot"
  
  
  # Construcción del biplot para Y (separado)
  
  BiplotY$nrows = I
  BiplotY$ncols = J*K
  BiplotY$dim = r1
  
  Y=BiplotY$Scaled_Data
  sct=sum(Y^2)
  
  scf = apply((Y^2), 1, sum)
  scc = apply((Y^2), 2, sum)
  
  a=t3pc$A
  rownames(a)=rownames(X)
  b=t(t3pc$H %*% kronecker( t(t3pc$C) , t(t3pc$B1) ))
  
  Esp=a%*%t(b)
  
  sce=rep(0,r1)
  cf=matrix(0,I,r1)
  cc=matrix(0,J*K,r1)
  for (i in 1:r1){
    Yesp=matrix(a[,i])%*%t(matrix(b[,i]))
    scfe = apply((Yesp^2), 1, sum)
    scce = apply((Yesp^2), 2, sum)
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
  scb = scb/(J*K)
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf
  
  rownames(b) = colnames(Y)
  colnames(a) = paste("Dim", 1:r1)
  
  BiplotY$RowCoordinates=a
  BiplotY$ColCoordinates=b
  
  
  rownames(cf) = rownames(Y)
  colnames(cf) = paste("Dim", 1:r1)
  BiplotY$RowContributions=cf
  
  rownames(cc) = colnames(Y)
  colnames(cc) = paste("Dim", 1:r1)
  BiplotY$ColContributions=cc
  
  BiplotY$Scale_Factor = scf
  BiplotY$ClusterType="us"
  BiplotY$Clusters = as.factor(matrix(1,nrow(BiplotY$RowContributions), 1))
  BiplotY$ClusterColors="blue"
  BiplotY$ClusterNames="Cluster 1"
  
  class(BiplotY)="ContinuousBiplot"
  
  t3pc$BiplotX=BiplotX
  t3pc$BiplotY=BiplotY
  t3pc$t3pcovr=t3pc
  return(t3pc)
}