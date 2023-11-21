Biplot.PLSRBIN <- function(plsr, BinBiplotType=1){
  X=plsr$X
  Y=plsr$Y

  I=dim(X)[1]
  J=dim(X)[2]
  K=dim(Y)[2]
  S=dim(plsr$XScores)[2]
  
  if (BinBiplotType==1){
  Biplot = list()
  Biplot$Title = " PLSR - Biplot"
  Biplot$Type = "PLSR-BR"
  Biplot$alpha=0
  Biplot$Dimension=S
  Biplot$Initial_Transformation=plsr$Initial_Transformation
  Biplot$ncols=J
  Biplot$nrows=I
  Biplot$dim=S
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  if (plsr$Initial_Transformation == "Within groups standardization")  Biplot$Deviations = plsr$Deviations
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile)[2, ]
  Biplot$P75 = apply(X, 2, quantile)[4, ]

  a=plsr$XScores
  b=plsr$XLoadings
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/I
  scb = scb/J
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf

  Biplot$RowCoordinates = a
  Biplot$ColCoordinates = b
  Biplot$Scale_Factor = scf
  apply(plsr$XScores^2, 2, sum)/sum(plsr$ScaledX^2)
  sum(plsr$ScaledX^2)
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
  class(Biplot)="ContinuousBiplot"

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
  Res$ColumnParameters=cbind(plsr$Intercepts, plsr$YLoadings)
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
  Res$Scale_Factor = scf
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
  Res$Tresholds = Res$ColumnParameters[, 1]/modelDif
  Res$Communalities = rowSums(Res$Loadings^2)
  Res$ColumnParameters[, 2:(S + 1)]=Res$ColumnParameters[, 2:(S + 1)]/scf
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
    
    
    a=plsr$XLoadings
    b=plsr$YLoadings
    sca = sum(a^2)
    scb = sum(b^2)
    sca = sca/I
    scb = scb/J
    scf = sqrt(sqrt(scb/sca))
    a = a * scf
    b = b/scf
    
    Biplot$RowCoordinates = a
    Biplot$ColCoordinates = b
    Biplot$Scale_Factor = scf
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
    class(Biplot)="ContinuousBiplot"
    
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
    Res$ColumnParameters=cbind(plsr$Intercepts, plsr$YLoadings)
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
    Res$Scale_Factor = scf
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
    Res$Tresholds = Res$ColumnParameters[, 1]/modelDif
    Res$Communalities = rowSums(Res$Loadings^2)
    Res$ColumnParameters[, 2:(S + 1)]=Res$ColumnParameters[, 2:(S + 1)]/scf
    Biplot$BinSupVarsBiplot=Res
    class(Biplot$BinSupVarsBiplot)="BinSupVarsBiplot"
  }
  return(Biplot)
}

