summary.Binary.Logistic.Biplot <- function(object, Normal=TRUE, latex=FALSE, kable=FALSE, digits=2, label="logbip", ...){
  
  dims=dim(object$ColumnParameters)[2]-1
  RR=cbind(object$Deviances, object$Dfs, object$pvalues, object$Nagelkerke, object$CoxSnell, object$MacFaden, object$PercentsCorrec*100, object$Sensitivity*100, object$Specificity*100)
  
  colnames(RR)=c("Deviance", "D.F", "P-val", "Nagelkerke", "Cox-Snell", "MacFaden", "% Correct", "Sensitivity", "Specificity")
  rownames(RR)=rownames(object$ColumnParameters)
  Total=c(object$DevianceTotal, object$TotalDf, object$p, object$TotNagelkerke, object$TotCoxSnell, object$TotMacFaden, object$TotalPercent*100, object$TotalSensitivity*100, object$TotalSpecificity*100)
  RR=rbind(RR,Total)
  
  LO=cbind(object$Tresholds, object$Loadings, object$Communalities)
  colnames(LO)=c("Thresholds", paste("Dim",1:dims,sep=""), "Communalities")
  rownames(LO)=rownames(object$ColumnParameters)
  
  eigen=apply(object$Loadings^2,2, sum)
  nvar=length(object$Loadings[,1])
  varia=cbind(eigen, eigen/nvar, cumsum(eigen)/nvar)
  rownames(varia)=paste("Dim",1:dims,sep="")
  colnames(varia)=c("Eigenvalue", "Percent", "Cummulative")
  
  if (Normal){
  print("BINARY LOGISTIC BIPLOT")
  print(paste("Type of Biplot : ", object$Type))
  print(paste("Initial Configuration : ", object$InitialConfig))
  print(paste("Method : ", object$Method))
  print(paste("Rotation : ", object$Rotation))
  print("-----------")
  
  
  print("COLUMN PARAMETERS")
  print(object$ColumnParameters)
  
  print("-----------")
  print("COLUMNS FIT")
  print(object$Deviances)

  print(RR)
  print("------------------------")
  print("Thresholds, Loadings and Communalities")
  print(LO)
  print("------------------------")
  print("Eigenvalues and Expained Variance")
  print(varia)
  cat("\n------------------------\n")
  cat("\n\n LATEX TABLES \n\n")
  }
  
  if (latex){
    cat("\n------------------------\n")
    cat("\n\n LATEX TABLES \n\n")
    print.xtable(xtable(RR, digits=digits, caption="Columns Fit measures"))
    print.xtable(xtable(LO, digits=digits, caption = "Thresholds, Loadings and Communalities"))
    print.xtable(xtable(varia, digits=digits, caption = "Eigenvalues and Expained Variance"))
  }
  
  if (kable){
    caption1=paste("\\label{tab:", label, "fit}", "Columns Fit measures", sep="")
    print(knitr::kable(RR, digits=digits, caption=caption1))
    caption2=paste("\\label{tab:", label, "loadings}", "Thresholds, Loadings and Communalities", sep="")
    print(knitr::kable(LO, digits=digits, caption = caption2))
    caption3=paste("\\label{tab:", label, "explained}", "Eigenvalues and Expained Variance", sep="")
    print(knitr::kable(varia, digits=digits, caption = caption3))
  }
 
}