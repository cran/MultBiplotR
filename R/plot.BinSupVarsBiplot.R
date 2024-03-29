plot.BinSupVarsBiplot <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                  ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, ColorVar="black", 
                                  WhatSupBinVars=NULL, CexScale=0.5,  ...){
  
  
  p=dim(x$ColumnParameters)[1]
  if (is.null(WhatSupBinVars)) WhatSupBinVars=rep(1,p)
  if (is.null(ColorVar)) ColorVar=rep("black", p)
  if (length(ColorVar)==1) ColorVar=rep(ColorVar, p)
  ColLabels=rownames(x$ColumnParameters)

  for (i in 1:p)
    if (WhatSupBinVars[i]){
    PlotBinaryVar(b0=x$ColumnParameters[i,1], bi1=x$ColumnParameters[i,F1+1], bi2=x$ColumnParameters[i,F2+1], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                  mode=mode, label=ColLabels[i], Color=ColorVar[i], CexScale=CexScale, ...)}
}
