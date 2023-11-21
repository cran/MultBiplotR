summary.TetraDualStatis <- function(object, ...) {
  
  cat("Call\n")
  print(object$call)
  cat("\n Correlations among occasions\n")
  print(object$RV)
  
  cat("\n Inter-Structure\n")
  Inercias=cbind(object$EigInter, object$InerInter, cumsum(object$InerInter))
  colnames(Inercias)=c("Eigenvalues", "Expl. Var.", "Cum.")
  print(Inercias)
  
  
  cat("\n \n \n Biplot for Statis------------\n")
  print(object$Biplot)
  
}