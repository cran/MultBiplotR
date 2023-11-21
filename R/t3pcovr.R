t3pcovr <- function( X , Y , I , J , K , L , r1=2 , r2=2 , r3=2 ,
                     conv=1e-6 , OriginalAlfa = 0.5, AlternativeLossF=1 ,
                     nRuns=100 , StartSeed=0)
{


  # The main input are matrices
  I=dim(X)[1]
  L=dim(X)[2]
  # You must provide a 3Way array convenienly prepared
  set.seed(StartSeed )
  #
  checkinput = 1

  if( r1 > min(I,J*K, L) )
  {
    cat(" ",fill=TRUE)
    cat("rank1 should be an integer between 1 and " , min(I,L) , fill=TRUE )
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if( r2 > min(J,I*K) )
  {
    cat(" ",fill=TRUE)
    cat("rank2 should be an integer between 1 and " , min(J,I*K) , fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if( r3 > min(K,I*J) )
  {
    cat(" ",fill=TRUE)
    cat("rank3 should be an integer between 1 and " , min(K,I*J) , fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  if ( (r1 > r2*r3) | (r2 > r1*r3) | (r3 > r1*r2) )
  {
    cat(" ",fill=TRUE)
    cat("None of the ranks can be larger than the products of the other two (e.g., rank1 > rank2*rank3 is not allowed)",fill=TRUE)
    cat(paste(r1,r2,r3))
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if ( (OriginalAlfa < 0) || (OriginalAlfa >= 1) )
  {
    cat(" ",fill=TRUE)
    cat("OriginalAlfa should be between 0 and 1 (but not 0 or 1)",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if ( (AlternativeLossF !=0) && (AlternativeLossF != 1) )
  {
    cat(" ",fill=TRUE)
    cat("AlternativeLossF should be 0 or 1",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }




  if ( checkinput == 1 )
  {
    ssq3D = sum(Y ^ 2)
    ssq2D = sum(X ^ 2)

    if (AlternativeLossF == 1)
    {
      Alfa = (OriginalAlfa * ssq3D) / ((OriginalAlfa * ssq3D) + ((1 - OriginalAlfa) * ssq2D))
    } else
    {
      Alfa = OriginalAlfa
    }
    BestA = matrix(0 , I , r1)
    BestB1 = matrix(0 , J , r2)
    BestC = matrix(0 , K , r3)
    BestH = array(0 , cbind(r1, r2, r3)) # CoreArray
    BestB2 = matrix(0 , L , r1)
    BestIter = -9999
    BestLoss = 999999999999

      FitValues = matrix(0 , 1 , nRuns + 1)
      nIterValues = matrix(0 , 1 , nRuns + 1)
      for (run in 1:nRuns + 1)
      {
        # initialize A, B1 and C (orthonormal)
        if (run  == 1)
        {
          # rational starts via eigendecomposition (gives orthonormal starting values)
          EIG = eigen(cbind(X, Y) %*% t(cbind(X, Y)))
          A = EIG$vectors[, 1:r1]
          rm(EIG)

          Z = permnew(Y , I , J , K)		# yields m x p x n array
          EIG = eigen(Z %*% t(Z))
          B1 = EIG$vectors[, 1:r2]
          rm(EIG)

          Z = permnew(Z , J , K , I)		# yields p x n x mrray
          EIG = eigen(Z %*% t(Z))
          C = EIG$vectors[, 1:r3]
          rm(EIG, Z)
        }
        else
        {
          # random start (orthonormal)
          A = orth(matrix(rnorm(I * r1 , 0 , 1) , I , r1))
          B1 = orth(matrix(rnorm(J * r2 , 0 , 1) , J , r2))
          C = orth(matrix(rnorm(K * r3 , 0 , 1) , K , r3))
        }
        # Calculate initial Core
        Z = permnew(t(A) %*% Y , r1 , J , K)
        Z = permnew(t(B1) %*% Z , r2 , K , r1)
        H = permnew(t(C) %*% Z , r3 , r1 , r2) # H (r1 x r2r3)
        rm(Z)

        # Calcule initial B2 (is in general not orthogonal !!! )
        B2 = t(ginv(t(A) %*% A) %*% t(A) %*% X)

        # Evaluate f
        Model3D = A %*% H %*% kronecker(t(C) , t(B1))
        Model2D = A %*% t(B2)
        f = (1 - Alfa) * sum((Y - Model3D) ^ 2)  +  Alfa * sum((X - Model2D) ^ 2)
        iter = 0
        fold = f + (2 * conv * f)
        V = cbind((sqrt(Alfa) * X) , (sqrt(1 - Alfa) * Y))
        while ((fold - f) > (f * conv))
        {
          iter = iter + 1
          fold = f

          U = t(cbind((sqrt(Alfa) * t(B2)), (
            sqrt(1 - Alfa) * H %*% kronecker(t(C) , t(B1))
          )))

          # update B1 (orthonormal)
          Z = permnew(Y , I , J , K)
          Z = permnew(Z , J , K , I)
          Z = permnew(t(C) %*% Z , r3 , I , J)
          Z = permnew(t(A) %*% Z , r1 , J , r3)			 # yields m x r3 x r1 array
          B1 = qr.Q(qr(Z %*% (t(Z) %*% B1)) , complete = FALSE)
          rm(Z)

          # update C (orthonormal)
          Z = permnew(t(A) %*% Y , r1 , J , K)
          Z = permnew(t(B1) %*% Z , r2 , K , r1)			 # yields p x r1 x r2 array
          C = qr.Q(qr(Z %*% (t(Z) %*% C)) , complete = FALSE)
          rm(Z)

          # Update H (Core)
          Z = permnew(t(A) %*% Y , r1 , J , K)
          Z = permnew(t(B1) %*% Z , r2 , K , r1)
          H = permnew(t(C) %*% Z , r3 , r1 , r2)
          rm(Z)

          #A=svd(X1 %*% K1)$u[,1:P]
          A = t(ginv(t(U) %*% U) %*% t(U) %*% t(V))
          A=orth(A)

          # Update B2 (not necessarily orthogonal !!)
          B2 = t(ginv(t(A) %*% A) %*% t(A) %*% X)


          # Evaluate f
          Model3D = A %*% H %*% kronecker(t(C) , t(B1))
          Model2D = A %*% t(B2)
          f = (1 - Alfa) * sum((Y - Model3D) ^ 2)  +  Alfa * sum((X - Model2D) ^ 2)
        }   #end of while-loop (alternating part of the algorithm)
        FitValues[run] = f
        nIterValues[run] = iter

        if (f < BestLoss)
        {
          BestA = A
          BestB1 = B1
          BestC = C
          BestB2 = B2
          BestH = H
          BestLoss = f
          BestIter = iter
        }
      } # end of for-loop (nRuns)

    W=MASS::ginv(t(X)%*%X)%*%t(X)%*%A
    # compute "intrinsic eigenvalues"
    La = BestH %*% t(BestH)
    J = permnew(BestH , r1 , r2 , r3)
    Lb = J %*% t(J)
    J = permnew(J , r2 , r3 , r1)
    Lc = J %*% t(J)

    # Compute BOF
    Model3D = BestA %*% BestH %*% kronecker(t(BestC) , t(BestB1))
    Model2D = BestA %*% t(BestB2)
    BOF3D = sum((Model3D - Y) ^ 2)
    BOF2D = sum((Model2D - X) ^ 2)

    FitPercentage = (((1 - Alfa) * (sum(BestH ^ 2) / ssq3D)) + (Alfa  * (sum(Model2D ^
                                                                               2) / ssq2D))) * 100
    FitPercentage3D = 100 * ((ssq3D - BOF3D) / ssq3D)
    FitPercentage2D = 100 * ((ssq2D - BOF2D) / ssq2D)

    out = list()
    out$Info = list()
    out$Info$nRows = I
    out$Info$nColumns3D = J
    out$Info$nSlices = K
    out$Info$nColumns2D = L
    out$Info$RankVector = cbind(r1 , r2 , r3)
    out$Info$TolPercentage = conv
    out$Info$OriginalAlfa = OriginalAlfa
    out$Info$Alfa = Alfa
    out$Info$AlternativeLossF = AlternativeLossF
    out$Info$nRuns = nRuns

    out$A = BestA
    out$B1 = BestB1
    out$C = BestC
    out$H = BestH
    out$B2 = BestB2
    out$LossWeighted = BestLoss
    out$LossUnweighted = BOF3D + BOF2D
    out$FitPercentage = FitPercentage
    out$FitPercentage3D = FitPercentage3D
    out$FitPercentage2D = FitPercentage2D
    out$nIter = BestIter
    out$FitValues = FitValues
    out$nIterValues = nIterValues
    out$La = La
    out$Lb = Lb
    out$Lc = Lc
    out$Fit3D = BOF3D
    out$Fit2D = BOF2D
    out$Fit = (Alfa * BOF3D) + ((1 - Alfa) * BOF2D)
    out$Fit3Dsize = BOF3D / (I * J * K)
    out$Fit2Dsize = BOF2D / (I * L)
    out$Fitsize = out$Fit / ((I * J * K) + (I * L))
    out$W=W
    out$Model3D=Model3D
    #out$CpuTime = cputime[1]
    #out$TimeSeconds = round(cputime[1] , 2)

    class(out) = "t3pcovr"
    return(out)
  }
}


