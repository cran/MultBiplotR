\name{Biplot.PLSRBIN}
\alias{Biplot.PLSRBIN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biplot for a PLSR model with binary responses
}
\description{
Builds a Biplot for a PLSR model with binary responses
}
\usage{
Biplot.PLSRBIN(plsr, BinBiplotType = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plsr}{
A PLSRBin object
}
  \item{BinBiplotType}{
The type of biplot:

   1:The biplot resulting from the fit,  for the binary responses.
   
   2: The biplot for the coefficients

}
}
\details{
Builds a Biplot for a PLSR model with binary responses. The result is a biplot for the matrix with the predictors (X) adding the binary responses as suplementary variables. There are two possible types, 1 for the biplot directly obtained in the fit ( the default) and 2 for the biplot obtaines after refitting the binary variables using Ridge Logistic Regression.
}
\value{
An object of class ContinuousBiplot  
}
\references{
Ugarte Fajardo, J., Bayona Andrade, O., Criollo Bonilla, R., Cevallos‐Cevallos, J., Mariduena‐Zavala, M., Ochoa Donoso, D., & Vicente Villardon, J. L. (2020). Early detection of black Sigatoka in banana leaves using hyperspectral images. Applications in plant sciences, 8(8), e11383.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
\donttest{
X=as.matrix(wine[,4:21])
Y=cbind(Factor2Binary(wine[,1])[,1], Factor2Binary(wine[,2])[,1])
rownames(Y)=wine[,3]
colnames(Y)=c("Year", "Origin")
pls=PLSRBin(Y,X, penalization=0.1, show=TRUE, S=2)
plsbip=Biplot.PLSRBIN(pls, BinBiplotType=1)
plsbip=AddCluster2Biplot(plsbip, ClusterType = "us", 
       Groups = wine$Group)
plot(plsbip, margin=0.05, mode="s", PlotClus = TRUE, 
    ModeSupBinVars = "s", ShowAxis = FALSE, 
    ColorSupBinVars = "blue",     CexInd=0.5, 
    ClustCenters = TRUE, LabelInd = FALSE, ShowBox = TRUE)
}
}
