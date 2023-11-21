BinaryPLSR <- function(Y, X, S=2, tolerance=0.00005,
                     maxiter=100, show=FALSE, penalization=0.1,
                     OptimMethod="CG", seed = 0){

  if (is.data.frame(X)) X=as.matrix(X)
  if (is.data.frame(Y)) X=as.matrix(Y)
  if (!CheckBinaryVector(Y)) stop("The response must be binary (0 or 1)")

  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering",
                              "Column centering", "Standardize columns", "Row centering",
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
 
  result=list()
  I1=dim(X)[1]
  J=dim(X)[2]

  I2=dim(Y)[1]
  K=1
  inames=rownames(X)
  ynames=colnames(Y)
  xnames=colnames(X)
  dimnames=paste("Comp.", 1:S)


  result$Method="PLSR for binary data"
  result$X=X
  result$Y=Y

  if (!(I1==I2)) stop('The number of rows of both matrices must be the same')
  else I=I1

  result$ScaledX=X
  result$ScaledY=Y
  result$tolerance=tolerance
  result$maxiter=maxiter
  result$penalization=penalization

  myfit=BinaryPLSFit(Y=Y, X=X, S=S, tolerance=tolerance, maxiter=maxiter, show=show, penalization=penalization, OptimMethod = OptimMethod, seed = seed)

  rownames(myfit$TT)=inames
  colnames(myfit$TT)=dimnames
  C=matrix(0, K, S)
  #rownames(myfit$B)=xnames
  #colnames(myfit$B)=ynames
  rownames(myfit$P)=xnames
  colnames(myfit$P)=dimnames

  result$XScores=myfit$TT
  result$XLoadings=myfit$P
  result$YScores=myfit$U
  result$YLoadings=myfit$Q
  rownames(result$YLoadings)=ynames
  colnames(result$YLoadings)=paste("Dim", 1:S)
  #result$Coefficients=myfit$B
  result$XStructure=cor(result$X,myfit$TT)
  result$BinaryFits=myfit$fit
  result$InterceptsY=myfit$q0
  result$InterceptsX=myfit$p0
  result$LinTerm=myfit$Linterm
  result$Expected=myfit$Expected
  result$Predictions=myfit$Predictions
  rownames(result$Predictions)=inames
  colnames(result$Predictions)=ynames
  result$PercentCorrect=myfit$PercentCorrect
  result$PercentCorrectCols=myfit$PercentCorrectCols
  result$OptimMethod=OptimMethod

  class(result)="BinaryPLSR"
  return(result)
}
