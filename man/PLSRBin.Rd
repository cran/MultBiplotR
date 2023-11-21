\name{PLSRBin}
\alias{PLSRBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial Least Squares Regression with several Binary Responses
}
\description{
Fits Partial Least Squares Regression with several Binary Responses
}

\usage{
PLSRBin(Y, X, S = 2, InitTransform = 5, grouping = NULL, 
tolerance = 5e-05, maxiter = 100, show = FALSE, penalization = 0.1, 
cte = TRUE, OptimMethod = "CG", Multiple = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The response
}
  \item{X}{
The matrix of independent variables
}
  \item{S}{
The Dimension of the solution
}
  \item{InitTransform}{
Initial transform for the X matrix
}
  \item{grouping}{
Grouping variable when the inial transformation is standardization within groups.
}
  \item{tolerance}{
Tolerance for convergence of the algorithm
}
  \item{maxiter}{
Maximum Number of iterations
}
  \item{show}{
Show the steps of the algorithm
}
  \item{penalization}{
Penalization for the Ridge Logistic Regression
}
  \item{cte}{
Should a constant be included in the model?
}
  \item{OptimMethod}{
Optimization methods from optim
}
  \item{Multiple}{
The responses are the indicators of a multinomial variable?
}
}
\details{
The function fits the PLSR method for the case when there is a set binary dependent variables, using logistic rather than linear fits to take into account the nature of responses. We term the method PLS-BLR (Partial Least Squares Binary Logistic Regression). This can be considered as a generalization of the NIPALS algorithm when the responses are all binary.
}
\value{
\item{Method }{Description of 'comp1'}
\item{X }{The predictors matrix}
\item{Y }{The responses matrix}
\item{Initial_Transformation }{Initial Transformation of the X matrix}
\item{ScaledX }{The scaled X matrix}
\item{tolerance }{Tolerance used in the algorithm}
\item{maxiter }{Maximum number of iterations used}
\item{penalization }{Ridge penalization}
\item{IncludeConst }{Is the constant included in the model?}
\item{XScores }{Scores of the X matrix, used later for the biplot}
\item{XLoadings }{Loadings of the X matrix}
\item{YScores }{Scores of the Y matrix}
\item{YLoadings }{Loadings of the Y matrix}
\item{Coefficients }{Regression coefficients}
\item{XStructure }{Correlations among the X variables and the PLS scores}
\item{Intercepts }{Intercepts for the Y loadings}
\item{LinTerm }{Linear terms for each response}
\item{Expected }{Expected probabilities for the responses}
\item{Predictions }{Binary predictions of the responses}
\item{PercentCorrect }{Global percent of correct predictions}
\item{PercentCorrectCols }{Percent of correct predictions for each column}
\item{Maxima }{Column with the maximum probability. Useful when the responses are the indicators of a multinomial variable}
}
\references{
Ugarte Fajardo, J., Bayona Andrade, O., Criollo Bonilla, R., Cevallos‐Cevallos, J., Mariduena‐Zavala, M., Ochoa Donoso, D., & Vicente Villardon, J. L. (2020). Early detection of black Sigatoka in banana leaves using hyperspectral images. Applications in plant sciences, 8(8), e11383.
}
\author{
José Luis Vicente Villardon
}

\examples{
\donttest{
X=as.matrix(wine[,4:21])
Y=cbind(Factor2Binary(wine[,1])[,1], Factor2Binary(wine[,2])[,1])
rownames(Y)=wine[,3]
colnames(Y)=c("Year", "Origin")
pls=PLSRBin(Y,X, penalization=0.1, show=TRUE, S=2)
}
}
