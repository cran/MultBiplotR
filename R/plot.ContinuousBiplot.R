# Plots a biplot for continuous data
plot.ContinuousBiplot <- function(x, A1 = 1, A2 = 2, ShowAxis = FALSE, margin = 0,
                                  PlotVars = TRUE, PlotInd = TRUE, WhatInds = NULL,
                                  WhatVars = NULL, LabelVars = TRUE, LabelInd = TRUE,
                                  IndLabels = NULL, VarLabels = NULL, mode = "a", CexInd
                                  = NULL, CexVar = NULL, ColorInd = NULL, ColorVar =
                                    NULL, LabelPos = 1, SmartLabels = FALSE, AbbreviateLabels=FALSE,
                                  MinQualityInds = 0, MinQualityVars = 0, dp = 0,
                                  PredPoints = 0, PlotAxis = FALSE, TypeScale =
                                    "Complete", ValuesScale = "Original", SizeQualInd =
                                    FALSE, SizeQualVars = FALSE, ColorQualInd = FALSE,
                                  ColorQualVars = FALSE, PchInd = NULL, PchVar = NULL,
                                  PlotClus = FALSE, TypeClus = "ch", ClustConf = 1,
                                  ClustLegend = FALSE, ClustLegendPos = "topright",
                                  ClustCenters = FALSE,  UseClusterColors = TRUE, CexClustCenters=1,
                                  PlotSupVars = TRUE, SupMode="a", ShowBox=FALSE, nticks=5, NonSelectedGray=FALSE, 
                                  PlotUnitCircle=TRUE, PlotContribFA=TRUE, AddArrow=FALSE, 
                                  ColorSupContVars="red", ColorSupBinVars="red", ColorSupOrdVars="red",
                                  ModeSupContVars="a", ModeSupBinVars="a", ModeSupOrdVars="a", 
                                  WhatSupBinVars=NULL, Title=NULL, Xlab=NULL, Ylab=NULL, add=FALSE,
                                  PlotTrajVars=FALSE, PlotTrajInds=FALSE, LabelTraj="end", Limits=NULL,
                                  PlotSupInds=FALSE, WhatSupInds=NULL, ColorSupInd="black", CexSupInd=0.8, PchSupInd=16,
                                  LabelSupInd=TRUE, PredSupPoints = 0, CexScale=0.5, ...){
  

  
  modes=c("p", "a", "b", "h", "ah", "s")
  if (is.numeric(mode)) 
    mode = modes[mode]
  TypeScales=c("Complete", "StdDev", "BoxPlot")
  if (is.numeric(TypeScale)) 
    TypeScale = TypeScales[TypeScale]
  ValuesScales=c("Original", "Transformed")
  if (is.numeric(ValuesScale)) 
    ValuesScale = ValuesScales[ValuesScale]
  
  # Obtaining coordinates and qualities for the representation
  A = x$RowCoordinates[, c(A1, A2)]
  B = x$ColCoordinates[, c(A1, A2)]
  n = dim(A)[1]
  p = dim(B)[1]
  
  if (is.null(IndLabels)) IndLabels=rownames(A)
  if (is.null(VarLabels)) VarLabels=rownames(B)
  
  if (AbbreviateLabels){
    IndLabels=abbreviate(IndLabels, minlength = 5L)
    VarLabels=abbreviate(VarLabels, minlength = 5L)
  }
  
  # Determining what rows to plot
  if (is.null(WhatInds)) 
    WhatInds = matrix(1, n, 1)
  else
    if (!CheckBinaryVector(WhatInds)) {
      AllRows = matrix(0, n, 1)
      AllRows[WhatInds]=1
      WhatInds=AllRows
    }
  WhatInds=as.logical(WhatInds)
  
  # Determining what columns to plot
  if (is.null(WhatVars)) 
    WhatVars = matrix(1, p, 1)
  else
    if (!CheckBinaryVector(WhatVars)){
      AllCols = matrix(0, p, 1)
      AllCols[WhatVars]=1
      WhatVars=AllCols
    }
  WhatVars=as.logical(WhatVars) 

  qlrcols = x$ColContributions[, A1] + x$ColContributions[, A2]
  qlrrows = x$RowContributions[, A1] + x$RowContributions[, A2]
  
  WhatInds = WhatInds & (qlrrows>MinQualityInds*100)
  WhatVars = WhatVars & (qlrcols>MinQualityVars*100)
  
  # Determining sizes and colors of the points
  if (is.null(CexInd)) 
    CexInd = rep(0.8, n)
  else if (length(CexInd == 1)) 
    CexInd = rep(CexInd, n)
  
  if (is.null(CexVar)) 
    CexVar = rep(0.8, p)
  else if (length(CexVar == 1)) 
    CexVar = rep(CexVar, p)
  
  if (is.null(PchInd)) 
    PchInd = rep(1,n)
  else if (length(PchInd == 1)) 
    PchInd = rep(PchInd, p)
  
  if (is.null(PchVar)) 
    PchVar = rep(16, p)
  else if (length(PchVar == 1)) 
    PchVar = rep(PchVar, p)
  
  if (is.null(ColorInd)) 
    if (is.null(x$ColorInd)) 
      ColorInd = rep("red",n)
  else ColorInd = x$ColorInd
  if (length(ColorInd)==1) ColorInd = rep(ColorInd,n)
  
  if (is.null(ColorVar)) 
    ColorVar = rep("black", p)
  if (length(ColorVar)==1) ColorVar = rep(ColorVar,p)
  if (SizeQualInd) 
    CexInd = cscale(qlrrows, rescale_pal())
  if (SizeQualVars) 
    CexVar = cscale(qlrcols, rescale_pal())
  
  if (NonSelectedGray){
    ColorInd[which(WhatInds==0)]="gray80"
    ColorVar[which(WhatVars==0)]="gray80"
    WhatInds = matrix(1, n, 1)
    WhatVars = matrix(1, p, 1)
  }
  
  
  if (ColorQualInd) 
    ColorInd = cscale(qlrrows, seq_gradient_pal("white", "red"))
  if (ColorQualVars) 
    ColorVar = cscale(qlrcols, seq_gradient_pal("white", "blue"))
  
  if (x$alpha == 2){ShowAxis=TRUE
  if(mode=="s") mode ="a"}
  
  if (ShowAxis) {
    xaxt = "s"
    yaxt = "s"
  } else {
    xaxt = "n"
    yaxt = "n"
  }
  
  if ((margin < 0) | (margin > 0.3)) 
    margin = 0
  
  
  if (PlotVars & PlotInd){
    P = rbind(A, B)
  }
  else
    if (PlotVars) P=B
  else P=A
  
  if (!LabelVars) VarLabels=rep(" ", p)
  
  if (is.null(Xlab)) Xlab=""
  if (is.null(Ylab)) Ylab=""
  
  if (!is.null(Title)){
    Xlab=paste("Dim", A1,"(",round(x$Inertia[A1], digits=1),"%)") 
    Ylab=paste("Dim", A2,"(",round(x$Inertia[A2], digits=1),"%)") 
  }
  
  if (is.null(Limits)){
    xmin = min(P[, 1])
    xmax = max(P[, 1])
    ymin = min(P[, 2])
    ymax = max(P[, 2])
  }
  else{
    xmin = Limits[1]
    xmax = Limits[2]
    ymin = Limits[3]
    ymax = Limits[4]
  }
    

  xrang=abs(xmax-xmin)
  yrang=abs(ymax-ymin)
  if (xmax <0 ) xmax=xmax*(-1)
  
  Xlims=c(xmin - (xmax - xmin) * margin, xmax + (xmax - xmin) * margin)
  
  Ylims=c(ymin - (ymax - ymin) * margin, ymax + (ymax - ymin) * margin)

#P = rbind(P, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin),c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
  
  if (!add)
  plot(P[, 1], P[, 2], cex = 0, asp = 1, xaxt = xaxt, yaxt = yaxt, xlab="", ylab="", bty="n",
       xlim=Xlims, ylim=Ylims, ...)
  #op=par(mai=c(0,0,0.5,0))
  #op=par(mar=c(1, 1, 1, 1) + 0.1)
  
  if (is.null(Title)){
  if (x$Type!="LogFreqBiplot")
    title(main = paste(x$Title,"(Dim", A1,"(",round(x$Inertia[A1], digits=1),"%)-",A2,"(",round(x$Inertia[A2], digits=1),"%))")  , omi = c(0, 0, 0, 0))
  else
    title(main = paste(x$Title,"(Dim", A1,"-",A2)  , omi = c(0, 0, 0, 0))}
  else{
    title(main = Title, omi = c(0, 0, 0, 0))
  }
 
  title(xlab=Xlab, ylab=Ylab, line=1)
  
  if (ShowBox) rect(xmin, ymin, xmax, ymax)
  
  if (ShowAxis) {abline(h = 0, v = 0, col = "gray60")
    box()}
  
  if (PlotClus) {
    ColorInd=PlotBiplotClusters(A, x$Clusters, TypeClus = TypeClus, ClusterColors = x$ClusterColors, 
                                ClusterNames=x$ClusterNames, centers = ClustCenters, ClustConf=ClustConf, 
                                CexClustCenters=CexClustCenters, Legend=ClustLegend, LegendPos= ClustLegendPos, ...)
    if (x$ClusterType=="gm"){
      ColorInd2=rgb((x$P %*% t(col2rgb(x$ClusterColors)))/255)
      if (UseClusterColors) ColorInd = ColorInd2
      PchInd=rep(16,n)
    }
  }
  

  
  
  if (PlotVars) {
    if ((x$Type == "FA") & (PlotContribFA)) AddArrow=TRUE
    if (grepl("s", mode)|grepl("i", mode))
      Scales = GetBiplotScales(x, nticks=nticks,  TypeScale = TypeScale, ValuesScale = ValuesScale)
    
  
    for (j in 1:p) 
      if (WhatVars[j]) 
      VarBiplot(B[j, 1], B[j, 2], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                label = VarLabels[j], mode = mode, CexPoint = CexVar[j], Color = ColorVar[j], 
                ticks = Scales$Ticks[[j]], ticklabels = Scales$Labels[[j]], ts = TypeScale, 
                PchPoint = PchVar[j], AddArrow=AddArrow, CexScale=CexScale, ...)
    
    for (idp in dp)
      if ((idp > 0) & (idp < (p + 1))) {
        g = B[idp, ]
        nn = (t(g) %*% g)
        scal <- (A %*% g)/nn[1, 1]
        Dscal <- diag(as.vector(scal))
        Fpr <- Dscal %*% matrix(rep(1, nrow(A)), ncol = 1) %*% t(g)
        nrFpr <- nrow(Fpr)
        dlines(A, Fpr, color=ColorVar[idp])
      }
    
    for (idp in PredPoints)
      if ((idp > 0) & (idp < (n + 1)))
        for (j in 1:p){
          g = B[j, ]
          nn = (t(g) %*% g)
          scal <- (A[idp,] %*% g)/nn[1, 1]
          Fpr <- scal %*% t(g)
          nrFpr <- nrow(Fpr)
          dlines(matrix(A[idp,],1,2) , Fpr, color=ColorVar[j])
        }
  }
  
  if (PlotTrajVars){
    points(B[, 1], B[, 2], type="l", lwd=2)
  }
  
  if (PlotTrajInds){
    if (LabelTraj=="end") lp=p else lp=1
    TrajInds=list()
    for (idp in 1:n){
      TrajInds[[idp]]=matrix(0,p,2)
      for (j in 1:p){
        g = B[j, ]
        nn = (t(g) %*% g)
        scal <- (A[idp,] %*% g)/nn[1, 1]
        Fpr <- scal %*% t(g)
        TrajInds[[idp]][j,]=Fpr
      }
      points(TrajInds[[idp]][,1], TrajInds[[idp]][,2], type="l", lwd=2, col=ColorInd[idp])
      text(TrajInds[[idp]][lp,1], TrajInds[[idp]][lp,2], col=ColorInd[idp], 
           label = IndLabels[idp], cex=0.7)
    }
  }
  
  if ((x$Type == "FA") & (PlotContribFA)){
    for (i in 1:10){
      Circle((i/10)/x$Scale_Factor, lty=3)
      #text(i/10, 0, labels=i/10, cex=0.5)
      text(0, (i/10)/x$Scale_Factor, labels=i/10, cex=0.5)
      #text(-1*i/10, 0, labels=i/10, cex=0.5)
      # text(0, -1*i/10, labels=i/10, cex=0.5)
    }
  }
  
  if (PlotSupVars) 
    plot.Supplementary.Variables(x, F1=A1, F2=A2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=SupMode, 
                                 TypeScale=TypeScale, ColorSupContVars=ColorSupContVars, ColorSupBinVars=ColorSupBinVars, 
                                 ColorSupOrdVars=ColorSupOrdVars, ModeSupContVars=ModeSupContVars, ModeSupBinVars=ModeSupBinVars, 
                                 ModeSupOrdVars=ModeSupOrdVars, WhatSupBinVars=WhatSupBinVars, CexScale=CexScale)
 
  
  if (PlotInd) 
    for (i in 1:n)
      if (WhatInds[i]){
        points(A[i, 1], A[i, 2], col = ColorInd[i], cex=CexInd[i], pch = PchInd[i], ...)
      }
  
  if (LabelInd)
    for (i in 1:n)
      if (WhatInds[i])
        if (SmartLabels) 
          textsmart(cbind(A[i, 1], A[i, 2]), CexPoints = CexInd[i], ColorPoints = ColorInd[i], ...)
  else text(A[i, 1], A[i, 2], IndLabels[i], cex = CexInd[i], col = ColorInd[i], pos = LabelPos, ...)
  
  if ((x$Type == "FA") & (PlotUnitCircle)){
    Circle(1/x$Scale_Factor)
  }
  
  
  if (PlotSupInds){
    nsi=dim(x$SupInds)[1]
    if (is.null(WhatSupInds)) 
      WhatSupInds = matrix(1, nsi, 1)
    else
      if (!CheckBinaryVector(PlotSupInds)) {
        AllRows = matrix(0, nsi, 1)
        AllRows[WhatSupInds]=1
        WhatSupInds=AllRows
      }
    WhatSupInds=as.logical(WhatSupInds)
    for (i in 1:nsi)
      if (WhatSupInds[i]){
        points(x$SupInd[i, 1], x$SupInd[i, 2], col = ColorSupInd, cex=CexSupInd, pch = PchSupInd, ...)
      }
    SupIndLabels=rownames(x$SupInd)
    if (LabelSupInd)
      for (i in 1:nsi)
      text(x$SupInd[i, 1], x$SupInd[i, 2], SupIndLabels[i], col = ColorSupInd, cex=CexSupInd, pos = LabelPos, ...)
    
    if (PlotVars) {
      for (idp in PredSupPoints)
        if ((idp > 0) & (idp < (nsi + 1)))
          for (j in 1:p){
            g = B[j, ]
            nn = (t(g) %*% g)
            scal <- (x$SupInds[idp,] %*% g)/nn[1, 1]
            Fpr <- scal %*% t(g)
            nrFpr <- nrow(Fpr)
            dlines(matrix(x$SupInds[idp,],1,2) , Fpr, color=ColorVar[j])
          }
    }
  }
  

  
    
  (par('usr'))
  # par(op)
}


