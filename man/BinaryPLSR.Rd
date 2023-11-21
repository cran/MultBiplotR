\name{BinaryPLSR}
\alias{BinaryPLSR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial Least Squares Regression with Binary Data
}
\description{
Fits Partial Least Squares Regression with Binary Data
}

\usage{
BinaryPLSR(Y, X, S = 2, tolerance = 5e-05, maxiter = 100, show = FALSE,
                   penalization = 0.1, OptimMethod = "CG", seed = 0)
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
  \item{OptimMethod}{
Optimization methods from optim
}
  \item{seed}{
Seed. By default is 0.
}
}
\details{
The function fits the PLSR method for the case when there are two sets of binary variables, using logistic rather than linear fits to take into account the nature of responses. We term the method BPLSR (Binary Partial Least Squares Regression). This can be considered as a generalization of the NIPALS algorithm when the data are all binary.
}
\value{
\item{Method }{Description of 'comp1'}
\item{X }{The predictors matrix}
\item{Y }{The responses matrix}
\item{ScaledX }{The scaled X matrix}
\item{tolerance }{Tolerance used in the algorithm}
\item{maxiter }{Maximum number of iterations used}
\item{penalization }{Ridge penalization}
\item{XScores }{Scores of the X matrix, used later for the biplot}
\item{XLoadings }{Loadings of the X matrix}
\item{YScores }{Scores of the Y matrix}
\item{YLoadings }{Loadings of the Y matrix}
%\item{Coefficients }{Regression coefficients}
\item{XStructure }{Correlations among the X variables and the PLS scores}
\item{InterceptsY}{Intercepts for the Y loadings}
\item{InterceptsX}{Intercepts for the Y loadings}
\item{LinTerm }{Linear terms for each response}
\item{Expected }{Expected probabilities for the responses}
\item{Predictions }{Binary predictions of the responses}
\item{PercentCorrect }{Global percent of correct predictions}
\item{PercentCorrectCols }{Percent of correct predictions for each column}
}
\references{
Ugarte Fajardo, J., Bayona Andrade, O., Criollo Bonilla, R., Cevallos‐Cevallos, J., Mariduena‐Zavala, M., Ochoa Donoso, D., & Vicente Villardon, J. L. (2020). Early detection of black Sigatoka in banana leaves using hyperspectral images. Applications in plant sciences, 8(8), e11383.

Vicente-Gonzalez, L., & Vicente-Villardon, J. L. (2022). Partial Least Squares Regression for Binary Responses and Its Associated Biplot Representation. Mathematics, 10(15), 2580.
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
