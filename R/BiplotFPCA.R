BiplotFPCA <- function(FPCA, X){

# desvs=apply(X, 2, sd)
# elim=which(desvs==0)
#  X=X[,-c(elim)]

  RowNames = rownames(X)
  VarNames = colnames(X)

  RowCoordinates=FPCA$scores
  rownames(RowCoordinates)=RowNames

  Biplot = list()
  Biplot$Title = "FPCA Biplot"
  Biplot$Type = "FPCA"
  Biplot$Non_Scaled_Data = X
  Biplot$alpha=1
  Biplot$Dimension=dim(RowCoordinates)[2]
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile, 0.25)
  Biplot$P75 = apply(X, 2, quantile, 0.75)
  Biplot$GMean = mean(X)
  Scaling="Column centering"
  Biplot$Initial_Transformation = Scaling
  Data = InitialTransform(X, transform = Scaling)
  X = Data$X
  rownames(X) = RowNames
  colnames(X) = VarNames
  Biplot$Scaled_Data = X
  Biplot$Structure=cor(X,RowCoordinates)

  b=t(solve(t(RowCoordinates)%*%RowCoordinates)%*%t(RowCoordinates)%*%X)

  n=dim(X)[1]
  p=dim(X)[2]
  Biplot$nrows = n
  Biplot$ncols = p
  Biplot$EigenValues = FPCA$varprop
  Biplot$Inertia = FPCA$varprop*100
  Biplot$CumInertia = cumsum(Biplot$Inertia)

  sca = sum(RowCoordinates^2)
  scb = sum(b^2)
  sca = sca/n
  scb = scb/p
  scf = sqrt(sqrt(scb/sca))
  RowCoordinates = RowCoordinates * scf
  b = b/scf

  Biplot$RowCoordinates=RowCoordinates
  Biplot$ColCoordinates=b
  Biplot$ColContributions=100*Biplot$Structure^2
  Biplot$RowContributions=matrix(1/Biplot$Dimension, n, Biplot$Dimension)
  Biplot$Scale_Factor = scf
  Biplot$ClusterType="us"
  Biplot$Clusters = as.factor(matrix(1,n, 1))
  Biplot$ClusterColors="blue"
  Biplot$ClusterNames="Cluster 1"
  class(Biplot) <- "ContinuousBiplot"
  return(Biplot)

}

