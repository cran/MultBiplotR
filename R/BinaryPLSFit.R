BinaryPLSFit <- function (Y, X, S = 2, tolerance = 5e-06, maxiter = 100, show = FALSE, 
                          penalization = 0.1, OptimMethod = "CG", seed = 0) 
{
  #Dimensiones X
  I1 = dim(X)[1]
  J = dim(X)[2]
  cat("Matrix X: ", I1, " rows and ", J, " columns \n")
  
  #Dimensiones Y
  I2 = dim(Y)[1]
  K = dim(Y)[2]
  cat("Matrix Y: ", I2, " rows and ", K, " columns \n")
  
  #Individuos
  I = I2
  
  # Iniciar matrices TT, P, Q y U
  TT = matrix(1, I, 1)
  P = matrix(0, J, S)
  Q = matrix(0, K, S)
  U=matrix(1,I,1)
  
  #Fijar semilla
  set.seed = seed
  
  #Calcular p0 y q0
  p0 = rnorm(J) 
  q0 = rnorm(K) 
  for (j in 1:J) p0[j] = RidgeBinaryLogistic(y = X[, j], matrix(1,I1, 1), penalization = 0)$beta
  for (j in 1:K) q0[j] = RidgeBinaryLogistic(y = Y[, j], matrix(1,I2, 1), penalization = 0)$beta
  P = p0
  Q = q0
  
  for (k in 1:S) {
    if (show) cat("\n Dimension:", k)
    tt = rnorm(I)
    TT = cbind(TT, tt)
    u = rnorm(I)
    U = cbind(U, u) 
    parP = rnorm(J)
    P = cbind(P, parP)
    parQ = rnorm(K)
    Q = cbind(Q, parQ)
    t = u
    error = 1
    iter = 0
    while ((error > tolerance) & (iter < maxiter)) {
      iter = iter + 1
      told = t
      resbipP <- optim(parP, fn = JLogBiplotRegBRec, gr = grLogBiplotRegBRec,
                       method = OptimMethod, X = X, A = TT, B = P, lambda = penalization)
      parP = resbipP$par
      P[, k + 1] = parP
      resbipT <- optim(t, fn = JLogBiplotRegARec, gr = grLogBiplotRegARec, 
                       method = OptimMethod, X = X, A = TT, B = P, lambda = penalization)
      t = resbipT$par
      t=scale(t)
      tt=t
      u=t
      U[,k+1]=u
      resbipQ <- optim(parQ, fn = JLogBiplotRegBRec, gr = grLogBiplotRegBRec, 
                       method = OptimMethod, X = Y, A = U, B = Q, lambda = penalization)
      parQ = resbipQ$par
      Q[, k + 1] = parQ
      resbipU <- optim(u, fn = JLogBiplotRegARec, gr = grLogBiplotRegARec, 
                       method = OptimMethod, X = Y, A = U, B = Q, lambda = penalization)
      u = resbipU$par
      U[, k + 1] = scale(u)
      t = u
      t=scale(t)
      TT[, k + 1] = t
      error = sum((told - t)^2)
      if (show) 
        cat("\n iteraction ", round(iter), ", rows ", round(J, 3), ", error ", round(error, 7), sep = "")
    }
    TT[, k+1] = tt
  }
  U = U[, -1]
  Q = Q[, -1]
  P = P[, -1]
  Lin1 = TT %*% t(cbind(q0, Q))
  Expected1 = exp(Lin1)/(1 + exp(Lin1))
  Pred = matrix(as.numeric(Expected1 > 0.5), nrow = I)
  Right = (Y == Pred)
  PercentCorrect = sum(Right)/(I * K)
  PercentCorrectCols = apply(Right, 2, sum)/I
  TT = TT[, -1]
  result = list(TT = TT, P = P, U = U, Q = Q, p0 = p0, q0 = q0, 
                Linterm = Lin1, Expected = Expected1, Predictions = Pred, 
                PercentCorrect = PercentCorrect, PercentCorrectCols = PercentCorrectCols)
  return(result)
}
