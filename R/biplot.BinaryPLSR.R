Biplot.BinaryPLSR <- function(plsr, BinBiplotType=1){
  X=plsr$X
  Y=plsr$Y

  I=dim(X)[1]
  J=dim(X)[2]
  K=dim(Y)[2]
  S=dim(plsr$XScores)[2]
  
  A = plsr$XScores
  B = cbind(plsr$InterceptsX,plsr$XLoadings)
  H=sigmoide(cbind(rep(1,I), A) %*% t(B))
  
  if (BinBiplotType==1){
  Biplot = list()
  Biplot$Title = "Binary PLSR - Biplot"
  Biplot$Type = "BPLSR"
  Biplot$alpha=0
  Biplot$Dimension=S
  Biplot$ncols=J
  Biplot$nrows=I
  Biplot$dim=S

  Biplot$Data=X
  Biplot$Penalization=plsr$penalization
  Biplot$Tolerance=plsr$tolerance
  Biplot$OptimMethod=plsr$OptimMethod
  Biplot$Biplot="Binary Logistic (Recursive Gradient Descent)"
  Biplot$Type= "Binary Logistic (Recursive Gradient Descent)"
  Biplot$RowCoordinates=A
  Biplot$ColumnParameters=B
  
  esp = cbind(rep(1,I), A) %*% t(B)
  pred = exp(esp)/(1 + exp(esp))
  pred2 = matrix(as.numeric(pred > 0.5), I, J)
  acier = matrix(as.numeric(round(X) == pred2), I, J)
  acierfil = 100*apply(acier,1,sum)/J
  aciercol = 100*apply(acier,2,sum)/I
  presences=apply(X, 2, sum)
  absences=I-presences
  sens = apply((acier==1) & (X==1), 2, sum)/presences
  spec = apply((acier==1) & (X==0), 2, sum)/absences
  totsens = sum((acier==1) & (X==1))/sum(presences)
  totspec = sum((acier==1) & (X==0))/sum(absences)
  gfit = (sum(sum(acier))/(I * J)) * 100
  esp0 = matrix(rep(1,I), I,1) %*% B[, 1]
  pred0 = exp(esp0)/(1 + exp(esp0))
  d1 = -2 * apply(X * log(pred0) + (1 - X) * log(1 - pred0),2,sum)
  d2 = -2 * apply(X * log(pred) + (1 - X) * log(1 - pred),2,sum)
  d = d1 - d2
  ps = matrix(0, J, 1)
  for (j in 1:J) ps[j] = 1 - pchisq(d[j], 1)
  
  Biplot$NullDeviances=d1
  Biplot$ModelDeviances=d2
  Biplot$Deviances=d
  Biplot$Dfs=rep(Biplot$dim, J)
  Biplot$pvalues=ps
  Biplot$CoxSnell=1-exp(-1*Biplot$Deviances/I)
  Biplot$Nagelkerke=Biplot$CoxSnell/(1-exp((Biplot$NullDeviances/(-2)))^(2/I))
  Biplot$MacFaden=1-(Biplot$ModelDeviances/Biplot$NullDeviances)
  Biplot$TotalPercent=gfit
  Biplot$ModelDevianceTotal=sum(Biplot$ModelDeviances)
  Biplot$NullDevianceTotal=sum(Biplot$NullDeviances)
  Biplot$DevianceTotal=sum(Biplot$Deviances)
  dd = sqrt(rowSums(cbind(1,Biplot$ColumnParameters[, 2:(Biplot$dim + 1)])^2))
  Biplot$Loadings = diag(1/dd) %*% Biplot$ColumnParameters[, 2:(Biplot$dim + 1)]
  rownames(Biplot$Loadings) = colnames(X)
  Biplot$Tresholds = Biplot$ColumnParameters[, 1]/d
  Biplot$Communalities = rowSums(Biplot$Loadings^2)
  nn=I*J
  Biplot$TotCoxSnell=1-exp(-1*Biplot$DevianceTotal/nn)
  Biplot$TotNagelkerke=Biplot$TotCoxSnell/(1-exp((Biplot$NullDevianceTotal/(-2)))^(2/nn))
  Biplot$TotMacFaden=1-(Biplot$ModelDevianceTotal/Biplot$NullDevianceTotal)
  Biplot$R2 = apply((X-H)^2,2, sum)/apply((X)^2,2, sum)
  Biplot$TotR2 = sum((X-H)^2) /sum((X)^2)
  pred= matrix(as.numeric(H>0.5),I , J)
  verdad = matrix(as.numeric(X==pred),I , J)
  Biplot$PercentsCorrec=apply(verdad, 2, sum)/I
  Biplot$TotalPercent=sum(verdad)/(I*J)
  Biplot$Sensitivity=sens
  Biplot$Specificity=spec
  Biplot$TotalSensitivity=totsens
  Biplot$TotalSpecificity=totspec
  Biplot$TotalDf = Biplot$dim*J
  Biplot$p=1-pchisq(Biplot$DevianceTotal, df = Biplot$TotalDf)
  Biplot$ClusterType="us"
  Biplot$Clusters = as.factor(matrix(1,I, 1))
  Biplot$ClusterColors="blue"
  Biplot$ClusterNames="ClusterTotal"
  class(Biplot) = "Binary.Logistic.Biplot"
  

  nulllinterm = matrix(1,I,1) %*% matrix(plsr$InterceptsY, nrow=1) 
  nullfitted = exp(nulllinterm)/(1 + exp(nulllinterm))
  nullDeviance = -2 * apply((Y * log(nullfitted) + (1 - Y) * log(1 - nullfitted)), 2,sum)
  modelDeviance = -2 * apply((Y * log( plsr$Expected) + (1 - Y) * log(1 -  plsr$Expected)), 2, sum)
  modelDif=(nullDeviance - modelDeviance)
  modeldf=S
  modelp=1-pchisq(modelDif, df =  modeldf)

  CoxSnell=1-exp(-1*modelDif/I)
  Nagelkerke=CoxSnell/(1-exp((nullDeviance/(-2)))^(2/I))
  MacFaden=1-(modelDeviance/nullDeviance)
  residuals=Y-plsr$Expected
  

  pred2 = matrix(as.numeric(plsr$Expected > 0.5), I, K)
  acier = matrix(as.numeric(round(Y) == pred2), I, K)
  presences=apply(Y, 2, sum)
  absences=I-presences
  sens = apply((acier==1) & (Y==1), 2, sum)/presences
  spec = apply((acier==1) & (Y==0), 2, sum)/absences
  totsens = sum((acier==1) & (Y==1))/sum(presences)
  totspec = sum((acier==1) & (Y==0))/sum(absences)
  
  
  SSRes=apply((residuals^2), 2,sum)
  SSTot=apply((Y^2), 2, sum)
  R2 = 1 - SSRes/SSTot
  
  Res=list()
  Res$Biplot="Binary Logistic (from PLS-BR)"
  Res$Type= "Binary Logistic (from PLS-BR)"
  Res$ColumnParameters=cbind(plsr$InterceptsY, plsr$YLoadings)
  rownames(Res$ColumnParameters)=rownames(plsr$YLoadings)
  #Res$ColumnParameters=Res$ColumnParameters
  Res$NullDeviances=nullDeviance
  Res$ModelDeviances=modelDeviance
  Res$ModelDevianceTotal=sum(Res$ModelDeviances)
  Res$NullDevianceTotal=sum(Res$NullDeviances)
  Res$Deviances=modelDif
  Res$Dfs=modeldf
  Res$pvalues=modelp
  Res$Bonferroni=(modelp * K)* ((modelp * K)<=1) + (((modelp * K)>1))
  Res$CoxSnell=CoxSnell
  Res$Nagelkerke=Nagelkerke
  Res$MacFaden=1-(Res$ModelDeviances/Res$NullDeviances)
  Res$R2=R2
  Prediction=matrix(as.numeric(plsr$Expected>0.5), I,K)
  Correct=matrix(as.numeric(Y==Prediction), I,K)
  Res$DevianceTotal=sum(Res$Deviances)
  nn=I*K
  Res$TotCoxSnell=1-exp(-1*Res$DevianceTotal/nn)
  Res$TotNagelkerke=Res$TotCoxSnell/(1-exp((Res$NullDevianceTotal/(-2)))^(2/nn))
  Res$TotMacFaden=1-(Res$ModelDevianceTotal/Res$NullDevianceTotal)
  Res$TotalDf=K*S
  SSResT=sum(residuals^2)
  SSTotT=sum(Y^2)
  Res$TotR2 = 1 - SSResT/SSTotT
  
  Res$PercentsCorrec=apply(Correct, 2, sum)/I
  Res$TotalPercent=sum(Correct)/(I*K)
  Res$Sensitivity=sens
  Res$Specificity=spec
  Res$TotalSensitivity=totsens
  Res$TotalSpecificity=totspec
  Res$p=1-pchisq(Res$DevianceTotal, df = Res$TotalDf)
  
  dd = sqrt(rowSums(cbind(1,Res$ColumnParameters[, 2:(S + 1)])^2))
  Res$Loadings = diag(1/dd) %*% Res$ColumnParameters[, 2:(S + 1)]
  rownames(Res$Loadings) = colnames(Y)
  Res$Tresholds = Res$ColumnParameters[, 1]/modelDif
  Res$Communalities = rowSums(Res$Loadings^2)
  #Res$ColumnParameters[, 2:(S + 1)]=Res$ColumnParameters[, 2:(S + 1)]/scf
  Biplot$BinSupVarsBiplot=Res
  class(Biplot$BinSupVarsBiplot)="BinSupVarsBiplot"
  }
  
  
  
  
  
  if (BinBiplotType==2){
    Biplot = list()
    Biplot$Title = "PLSR - Biplot of the coefficients"
    Biplot$Type = "PLSR-BR Coef"
    Biplot$alpha=0
    Biplot$Dimension=S
    Biplot$Initial_Transformation="None"
    Biplot$ncols=K
    Biplot$nrows=J
    Biplot$dim=S
    
    
    a=cbind(plsr$InterceptsX,plsr$XLoadings)
    b=cbind(plsr$InterceptsY,plsr$YLoadings)
    
    Biplot$RowCoordinates = a
    Biplot$ColCoordinates = b
    apply(plsr$XScores^2, 2, sum)/sum(plsr$ScaledX^2)
    
    Cont=CalculateContributions(plsr$ScaledX, plsr$XScores, plsr$XLoadings )
    Biplot$EigenValues=apply(plsr$XScores^2, 2, sum)
    Biplot$Inertia=Cont$Fit*100
    Biplot$CumInertia=cumsum(Biplot$Inertia)
    Biplot$RowContributions=Cont$RowContributions
    Biplot$Structure=Cont$Structure
    Biplot$ColContributions=Cont$ColContributions
    Biplot$CumRowContributions=Cont$RowContributions
    Biplot$CumColContributions=Cont$CumColContributions
    Biplot$SupInds=b
    class(Biplot)="Binary.Logistic.Biplot"
    
    nulllinterm = matrix(1,I,1) %*% matrix(plsr$Intercepts, nrow=1) 
    nullfitted = exp(nulllinterm)/(1 + exp(nulllinterm))
    nullDeviance = -2 * apply((Y * log(nullfitted) + (1 - Y) * log(1 - nullfitted)), 2,sum)
    modelDeviance = -2 * apply((Y * log( plsr$Expected) + (1 - Y) * log(1 -  plsr$Expected)), 2, sum)
    modelDif=(nullDeviance - modelDeviance)
    modeldf=S
    modelp=1-pchisq(modelDif, df =  modeldf)
    
    CoxSnell=1-exp(-1*modelDif/I)
    Nagelkerke=CoxSnell/(1-exp((nullDeviance/(-2)))^(2/I))
    MacFaden=1-(modelDeviance/nullDeviance)
    residuals=Y-plsr$Expected
    
    
    pred2 = matrix(as.numeric(plsr$Expected > 0.5), I, K)
    acier = matrix(as.numeric(round(Y) == pred2), I, K)
    presences=apply(Y, 2, sum)
    absences=I-presences
    sens = apply((acier==1) & (Y==1), 2, sum)/presences
    spec = apply((acier==1) & (Y==0), 2, sum)/absences
    totsens = sum((acier==1) & (Y==1))/sum(presences)
    totspec = sum((acier==1) & (Y==0))/sum(absences)
    
    
    SSRes=apply((residuals^2), 2,sum)
    SSTot=apply((Y^2), 2, sum)
    R2 = 1 - SSRes/SSTot
    
    Res=list()
    Res$Biplot="Binary Logistic (from PLS-BR)"
    Res$Type= "Binary Logistic (from PLS-BR)"
    Res$ColumnParameters=cbind(plsr$InterceptsY, plsr$YLoadings)
    rownames(Res$ColumnParameters)=rownames(plsr$YLoadings)
    #Res$ColumnParameters=Res$ColumnParameters
    Res$NullDeviances=nullDeviance
    Res$ModelDeviances=modelDeviance
    Res$ModelDevianceTotal=sum(Res$ModelDeviances)
    Res$NullDevianceTotal=sum(Res$NullDeviances)
    Res$Deviances=modelDif
    Res$Dfs=modeldf
    Res$pvalues=modelp
    Res$Bonferroni=(modelp * K)* ((modelp * K)<=1) + (((modelp * K)>1))
    Res$CoxSnell=CoxSnell
    Res$Nagelkerke=Nagelkerke
    Res$MacFaden=1-(Res$ModelDeviances/Res$NullDeviances)
    Res$R2=R2
    Prediction=matrix(as.numeric(plsr$Expected>0.5), I,K)
    Correct=matrix(as.numeric(Y==Prediction), I,K)
    Res$DevianceTotal=sum(Res$Deviances)
    nn=I*K
    Res$TotCoxSnell=1-exp(-1*Res$DevianceTotal/nn)
    Res$TotNagelkerke=Res$TotCoxSnell/(1-exp((Res$NullDevianceTotal/(-2)))^(2/nn))
    Res$TotMacFaden=1-(Res$ModelDevianceTotal/Res$NullDevianceTotal)
    Res$TotalDf=K*S
    SSResT=sum(residuals^2)
    SSTotT=sum(Y^2)
    Res$TotR2 = 1 - SSResT/SSTotT
    
    Res$PercentsCorrec=apply(Correct, 2, sum)/I
    Res$TotalPercent=sum(Correct)/(I*K)
    Res$Sensitivity=sens
    Res$Specificity=spec
    Res$TotalSensitivity=totsens
    Res$TotalSpecificity=totspec
    Res$p=1-pchisq(Res$DevianceTotal, df = Res$TotalDf)
    
    dd = sqrt(rowSums(cbind(1,Res$ColumnParameters[, 2:(S + 1)])^2))
    Res$Loadings = diag(1/dd) %*% Res$ColumnParameters[, 2:(S + 1)]
    rownames(Res$Loadings) = colnames(Y)
    Res$Tresholds = Res$ColumnParameters[, 1]/modelDif
    Res$Communalities = rowSums(Res$Loadings^2)
    #Res$ColumnParameters[, 2:(S + 1)]=Res$ColumnParameters[, 2:(S + 1)]/scf
    Biplot$BinSupVarsBiplot=Res
    class(Biplot$BinSupVarsBiplot)="BinSupVarsBiplot"
  }
  return(Biplot)
}

