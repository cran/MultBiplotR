summary.ContinuousBiplot <- function(object, latex=FALSE, ...) {
  if (object$Type=="PCA")
    cat(" ###### Biplot for Principal Components Analysis #######\n\n")
  if (object$Type=="FA")
    cat(" ###### Biplot for Factor Analysis #######\n\n")
  if (object$Type=="HJ")
    cat(" ###### HJ - Biplot #######\n\n")
  if (object$Type=="PLSR")
    cat(" ###### Biplot for PLS Regression #######\n\n")
    if (object$Type=="PLSR_BR")
    cat(" ###### Biplot for PLS Binary Regression #######\n\n")
  
  
  cat("Call\n")
  print(object$call)
  cat("Type of coordinates:\n")  
  if (object$alpha==2)
    tipo="Principal Normalization (Baricentric Scaling)"
  if (object$alpha==1)
    tipo="Row Principal Normalization (RMP-Biplot)"
  if (object$alpha==0)
    tipo="Column Principal Normalization (CMP-Biplot)"
  if (object$alpha==0.5)
    tipo="Symmetrical Normalization (SQRT - Biplot)"
  if ((object$alpha>0) & (object$alpha<1) & (object$alpha != 0.5))
    tipo=paste("Custom Normalization (Biplot con \alpha = ",gamma,")")
  cat("Transformation of the raw data:\n")
  print(object$Initial_Transformation)
  cat("Type of Biplot\n")
  print(object$Type)
  cat("\n Eigenvalues & Explained Variance (Inertia)\n")
  
  pp=cbind(object$EigenValues[1:object$Dimension], object$Inertia[1:object$Dimension], object$CumInertia[1:object$Dimension])
  
  colnames(pp)=c("Eigenvalue", "Exp. Var", "Cummulative")
  print(pp)
  
  if (object$Type=="FA"){
    cat(" Factor Analysis Loadings\n\n")
    print(object$ColCoordinates)
  }
  
  cat("\n\n RELATIVE CONTRIBUTIONS OF THE FACTOR TO THE ELEMENT\n")
  cat("\n Row Contributions \n")
  print(round(object$RowContributions, digits=2))
  cat("\n Column Contributions \n")
  print(round(object$ColContributions, digits=2))
  

  cat("\n\n\n Qualities of representation of the rows (Cummulative contributions) \n")
  if (is.null(object$CumRowContributions))
    print(round(t(apply(object$RowContributions,1, cumsum)), digits=2))
  else
    print(round(object$CumRowContributions, digits=2))
  
  cat("\n\n\n Qualities of representation of the columns (Cummulative contributions) \n")
  if (is.null(object$CumColContributions))
    print(round(t(apply(object$ColContributions,1, cumsum)), digits=2))
  else
    print(round(object$CumColContributions, digits=2))
  
  cat("\n\n\n LATEX TABLES \n\n")
  
  if (latex){
    print(xtable(pp, caption="Explained Variance"))
    print(xtable(round(object$RowContributions, digits=2), caption="Row Contributions Factor to element"))
    print(xtable(round(object$ColContributions, digits=2), caption="Column Contributions Factor to element"))
    
    if (is.null(object$CumRowContributions))
      print(xtable(round(t(apply(object$RowContributions,1, cumsum)), digits=2), caption="Qualities of the rows"))
    else
      print(xtable(round(object$CumRowContributions, digits=2), caption="Qualities of the rows"))
    
    if (is.null(object$CumColContributions))
      print(xtable(round(t(apply(object$ColContributions,1, cumsum)), digits=2), caption="Qualities of the columns"))
    else
      print(xtable(round(object$CumColContributions, digits=2), caption="Qualities of the columns"))
    
  }
  
  if (!is.null(object$BinSupVarsBiplot)){
    cat("\n\n\ SUPPLEMENTARY BINARY VARIABLES\n")
    summary(object$BinSupVarsBiplot, latex=latex, ...)}
}
