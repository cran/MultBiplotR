plot.Supplementary.Variables <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                         ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, 
                                         ColorSupContVars="black", ColorSupBinVars="darkblue", ColorSupOrdVars="blue", 
                                         ModeSupContVars="s", ModeSupBinVars="s", 
                                         ModeSupOrdVars="s", WhatSupBinVars=NULL, CexScale=0.5, ...){
  

   if (!is.null(x$ContSupVarsBiplot)){
     if (is.null(ModeSupContVars)) ModeSupContVars=mode
    plot(x$ContSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=ModeSupContVars, TypeScale=TypeScale, ColorVar=ColorSupContVars, CexScale=CexScale) 
  }
  if (!is.null(x$BinSupVarsBiplot))
    plot(x$BinSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, 
         ymin = ymin, ymax = ymax, mode=ModeSupBinVars, 
         ColorVar=ColorSupBinVars, WhatSupBinVars=WhatSupBinVars, CexScale=CexScale)

  if (!is.null(x$OrdSupVarsBiplot))
    plot(x$OrdSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=ModeSupOrdVars, ColorVar=ColorSupOrdVars, CexScale=CexScale) 
  
}
