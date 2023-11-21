Sim.Biplot.t3pcovr <- function(sol, TypeBiplot=1){
  # Biplot Types
  # 1.- Prediction for X and Y
  # 2.- Regression coefficients
  # 3.- Interpolation for X and prediction for Y
  
  
  if (TypeBiplot==1){
    B=t(sol$H %*% kronecker( t(sol$C) , t(sol$B1) ))
    A=sol$A
    Biplot=sol$BiplotY
    sca = sum(A^2)
    scb = sum(B^2)
    sca = sca/dim(A)[1]
    scb = scb/dim(B)[1]
    scf = sqrt(sqrt(scb/sca))
    A = A * scf
    B = B/scf
    Biplot$RowCoordinates=A
    Biplot$ColCoordinates=B
    CorrXCP = cor(Biplot$Scaled_Data, A, method = "pearson")
    
    Biplot$nrows = dim(A)[1]
    Biplot$ncols = dim(B)[1]
    # Inertia
    Biplot$dim = dim(A)[2]
    
  }

  Yestimada =A %*% t(B)
  
  100*apply(Yestimada^2, 1, sum)/apply(Biplot$Scaled_Data^2, 1, sum)
  100*apply(Yestimada^2, 2, sum)/apply(Biplot$Scaled_Data^2, 2, sum)
  
  resid=Biplot$Scaled_Data-Yestimada
  sum(resid^2)
  sum(Yestimada^2)/sum(Biplot$Scaled_Data^2) 
  

  class(Biplot)="ContinuousBiplot"
  
  return(Biplot)
}