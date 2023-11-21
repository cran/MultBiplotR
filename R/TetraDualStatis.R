TetraDualStatis <- function(X, dimens=2, SameInd=FALSE, RotVarimax=FALSE, OptimMethod="L-BFGS-B",
                            penalization=0.01) {
  mycall=match.call()
  
  ng = length(X) #Number of groups
  StudyNames=names(X)
  nri = matrix(0, ng, 1) #Number of rows in each group
  nci = matrix(0, ng, 1) #Number of cols in each group
  for (i in 1:ng) {
    nri[i] = dim(X[[i]])[1]
    nci[i] = dim(X[[i]])[2]
  }
  nc = nci[1]
  if (sum(nci == nc) < ng) 
    stop("The number of columns (variables) must be the same in all ocassions")
  
  nr=sum(nri)

  #  Extracting the names of the occasions
  OccNames=names(X)
  if (is.null(OccNames)) {
    for (i in 1:ng) OccNames= c(OccNames, paste("Occasion_",i,sep=""))
  }
  
  # Initial transformation of data and calculation of statistics for the biplot
  
  StatisRes=list()
  StatisRes$Title="Dual STATIS Tetrachoric"
  StatisRes$Type="Tetra Dual STATIS"
  StatisRes$NTables=ng
  StatisRes$NCols=nc
  StatisRes$NVars=nc
  
  if (SameInd)
  {StatisRes$NRows=nr}
  else
  {StatisRes$NRows=nri}
  
  StatisRes$SameVar=TRUE
  
  
  StatisRes$VarLabels=colnames(X[[1]])
  StatisRes$TableLabels=names(X)

  #Calculation of the objects
  Ct = list()
  Taut=list()
  Tau=matrix(0, nc, ng)
  for (i in 1:ng) {
    Tet=tetrachoric(X[[i]])
    Ct[[i]] = Tet$rho
    Taut[[i]]=Tet$tau
    Tau[,i]=Tet$tau
  }
  
  # Calculo de los productos escalares
  S = matrix(0, ng, ng)
  for (i in 1:ng) for (j in 1:ng) {
    S[i, j] = sum(diag(Ct[[i]] %*% Ct[[j]]))
    S[j, i] = S[i, j]
  }
  rownames(S)=StudyNames
  colnames(S)=StudyNames
  
  #RV: correlations among the occasions
  RV = solve(diag(sqrt(diag(S)))) %*% S %*% solve(diag(sqrt(diag(S))))
  rownames(RV)=StudyNames
  colnames(RV)=StudyNames
  
  Inter = svd(RV)
  Weights = abs(Inter$u[, 1])/sum(abs(Inter$u[, 1]))
  
  # This information will be common to all the techniques, then a single
  # function to plot is necassary
  StatisRes$Data = X
  StatisRes$Objects = Ct
  StatisRes$Taus=Taut
  StatisRes$ScalarProducts = S
  StatisRes$RV = RV
  StatisRes$EigInter = Inter$d
  names(StatisRes$EigInter)=paste("Dim",1:length(StudyNames))
  
  StatisRes$InerInter=(Inter$d/sum(Inter$d))*100
  names(StatisRes$InerInter)=paste("Dim",1:length(StudyNames))
  StatisRes$InterStructure = -1 * Inter$u %*% diag(sqrt(Inter$d))
  rownames(StatisRes$InterStructure)=StudyNames
  colnames(StatisRes$InterStructure)=paste("Dim", 1:ng)
  
  # Weigths for the Compromise
  StatisRes$Weights = Weights
  
  # Compromiso para las correlaciones y los umbrales
  C=matrix(0,nc,nc)
  Taus=matrix(0, 1, nci)
  for (i in 1:ng){
    C=C + StatisRes$Weights[i]*Ct[[i]]
    Tau=Taus+StatisRes$Weights[i]*Taut[[i]]
    }
  StatisRes$C <- C
  colnames(Taus)=names(Taut[[1]])
  StatisRes$Tau <- Taus
  
  # Euclidean Configuration of the compromise
  Intra = svd(C)
  r=sum(Intra$d>0.00001)
  
  P=  -1 * Intra$u[,1:dimens] %*% diag(sqrt(Intra$d[1:dimens]))
 
  if (RotVarimax) {
    BB=varimax(P)
    P= P %*% BB$rotmat
    eigenval=apply(BB$loadings^2, 2, sum)
    Intra$d[1:dimens]=eigenval
  } 
  
  # El el primer biplot (el factorial) las coordenadas estan en 
  
  BiplotStatis=list()
  BiplotStatis$call <- mycall
  BiplotStatis$Type=paste("DualStatis", "Binary")
  BiplotStatis$Title="Tetrachoric Dual-Statis"
  BiplotStatis$nrows=sum(nri)
  BiplotStatis$ncols=nc
  BiplotStatis$EigenValues=Intra$d[1:r]
  BiplotStatis$Inertia=(Intra$d[1:r]/sum(diag(Intra$d[1:r])))*100
  BiplotStatis$CumInertia = cumsum(BiplotStatis$Inertia)
  BiplotStatis$alpha=0
  BiplotStatis$Dimension=dimens

  
  sf=apply((Intra$u[,1:r] %*% diag(sqrt(Intra$d[1:r])))^2,1,sum)
  BiplotStatis$ColContributions=(diag(1/sf) %*% P^2)*100
  rownames(BiplotStatis$ColContributions)=colnames(X[[1]])
  colnames(BiplotStatis$ColContributions)=paste("Dim", 1:dimens)
  BiplotStatis$Structure=P
  rownames(BiplotStatis$Structure)=colnames(X[[1]])
  colnames(BiplotStatis$Structure)=paste("Dim", 1:dimens)
  
  escala=1/sqrt(1+apply(P^2, 1, sum))
  b0=Taus*escala
  B=diag(escala)%*%P
  B=cbind(t(b0),B)
  
  CoorRows=list()
  RowCoordinates=NULL
  rowcolors=NULL
  for (i in 1:ng){
    parA=runif(nri[i]*dimens)
    XX=as.matrix(X[[i]])
    A=matrix(parA,nri[i],dimens)
    resbipA <- optim(parA, fn=JLogBiplotRegA, gr=grLogBiplotRegA, 
                      method=OptimMethod, X=XX, B=B, lambda=penalization)
    parA=resbipA$par
    A=matrix(parA,nri[i],dimens)
    CoorRows[[i]]=A
    RowCoordinates=rbind(RowCoordinates, A)
    rowcolors=c(rowcolors, rep(i, nri[i]))
  }
  
  
  
  
  # Reescalado de las coordenads para el biplot
  
  sca = sum(RowCoordinates^2)
  scb = sum(P^2)
  sca = sca/sum(nri)
  scb = scb/nc
  scf = sqrt(sqrt(scb/sca))
  RowCoordinates = RowCoordinates * scf
  P = P/scf
  
  BiplotStatis$RowCoordinates=RowCoordinates
  BiplotStatis$RowColors=rowcolors
  
  BiplotStatis$ColCoordinates = P
  rownames(BiplotStatis$ColCoordinates)=colnames(X[[1]])
  colnames(BiplotStatis$ColCoordinates)=paste("Dim", 1:dimens)
  
  BiplotStatis$RowContributions=matrix(1, sum(nri), dimens)
  trajvar=list()
  for (j in 1:ng){
    trajvar[[j]]=-1*Ct[[j]] %*% Intra$u[,1:dimens]  %*% diag(sqrt(1/Intra$d[1:dimens]))
    if (RotVarimax) {
      trajvar[[j]]= trajvar[[j]] %*% BB$rotmat
    } 
  }
  # Coordinates of the variables on the biplot
  
  
  
  BiplotStatis$TrajVar=trajvar
  
  class(BiplotStatis) ="ContinuousBiplot"
  StatisRes$Biplot=BiplotStatis
  class(StatisRes) = "TetraDualStatis"
  return(StatisRes)
}
