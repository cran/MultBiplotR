PolyOrdinalLogBiplot <- function (X, dimension = 3, method="principal", rotate="varimax", 
                                  RescaleCoordinates=TRUE, ...) {
# X must be a mtrix with integer numbers with the categories of the ordinal variables
# Method can be c("principal", "fa")
  if (is.data.frame(X)) 
    X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  # Setting the properties of data
  if (is.null(rownames(X))) 
    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "I")
  RowNames = rownames(X)
  if (is.null(colnames(X))) 
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "V")
  VarNames = colnames(X)
  
  DimNames = paste("Factor", 1:dimension, sep="_")
  
  # Polychoric correlation
  R=polychoric(X)
  Tau=R$tau
  
  if (method=="principal")
    res.fa=principal(X, nfactors=dimension, rotate = rotate, covar=FALSE,  cor="poly", ...)
  else
    res.fa=fa(X, nfactors=dimension, rotate = rotate, covar=FALSE,  cor="poly", ...)
  
  # res.fa=fa(X, nfactors=3, rotate = rotate, covar=FALSE,  cor="poly")
  
  a=as.matrix(res.fa$scores)
  b=as.matrix(res.fa$loadings)
  
  Biplot = list()
  Biplot$Title = "Factor Analysis Biplot (Ordinal Data)"
  Biplot$Type = "FA" 
  Biplot$call <- match.call()
  Biplot$alpha=0
  Biplot$Dimension=dimension
  Biplot$Non_Scaled_Data = X
  Biplot$Means = rep(0,p)
  Biplot$Medians = rep(0,p)
  Biplot$Deviations = rep(1,p)
  Biplot$Minima = rep(-3,p)
  Biplot$Maxima = rep(3,p)
  Biplot$P25 = rep(-0.67,p)
  Biplot$P75 = rep(0.67,p)
  Biplot$GMean = 0
  Biplot$Initial_Transformation = 1
  Biplot$Scaled_Data = X
  Biplot$R = R
  
  EV=apply(b^2,2,sum)
  Inertia = round((EV/p) * 100, digits = 3)
  CumInertia = cumsum(Inertia)
  Biplot$Communalities=apply(b^2,1,sum)
  Biplot$Uniqueness=1-Biplot$Communalities
  rownames(a) <- RowNames
  colnames(a) <- DimNames
  rownames(b) <- VarNames
  colnames(b) <- DimNames
  
  # Es necesario revisar las contribuciones de las filas. No se cual es la 
  # acotación en el espacio completo. Lo que está claro es que no es la suma
  # de cuadrados como en el caso continuo
  
  
  sf = apply((X^2), 1, sum)
  cf=matrix(0,n,dimension)
  for (k in 1:dimension)
    cf[,k]= round((a[,k]*sqrt(EV[k]))^2/ sf*100, digits = 2)
  rownames(cf) = RowNames
  colnames(cf) = DimNames
  cfacum = t(apply(cf, 1, cumsum))
  # Relative contributions of the rows 
  
  cc=round((b^2)*100, digits = 2)
  rownames(cc) = VarNames
  colnames(cc) = DimNames
  ccacum = t(apply(cc, 1, cumsum))
  
  
  scf = 1
  
  if (RescaleCoordinates){
    sca = sum(a^2)
    scb = sum(b^2)
    sca = sca/n
    scb = scb/p
    scf = sqrt(sqrt(scb/sca))
    a = a * scf
    b = b/scf
  }
  
  
  Biplot$nrows = n
  Biplot$ncols = p
  Biplot$dim = dimension
  Biplot$EigenValues = EV
  Biplot$Inertia = Inertia
  
  Biplot$CumInertia = CumInertia
  Biplot$Structure=unclass(res.fa$Structure)
  
  Biplot$RowCoordinates <- a
  Biplot$ColCoordinates <- b

  # Contributions

  Biplot$RowContributions <- cf
  Biplot$ColContributions <- unclass(cc)
  
  Biplot$Scale_Factor = scf
  Biplot$ClusterType="us"
  Biplot$Clusters = as.factor(matrix(1,n, 1))
  Biplot$ClusterColors="blue"
  Biplot$ClusterNames="Cluster 1"
  class(Biplot) <- "ContinuousBiplot"
  return(Biplot)

}
