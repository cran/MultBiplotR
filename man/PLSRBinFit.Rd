\name{PLSRBinFit}
\alias{PLSRBinFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PLS binary regression.
}
\description{
Fits PLS binary regression.
}
\usage{
PLSRBinFit(Y, X, S = 2, tolerance = 5e-06, maxiter = 100, 
show = FALSE, penalization = 0.1, cte = TRUE, OptimMethod = "CG")
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
  \item{cte}{
Should a constant be included in the model?
}
  \item{OptimMethod}{
Optimization methods from optim
}
}
\details{
Fits PLS binary regression. It is used for a higher level function.
}
\value{
The PLS fit used by the PLSRBin function.
}
\references{
Ugarte Fajardo, J., Bayona Andrade, O., Criollo Bonilla, R., Cevallos‐Cevallos, J., Mariduena‐Zavala, M., Ochoa Donoso, D., & Vicente Villardon, J. L. (2020). Early detection of black Sigatoka in banana leaves using hyperspectral images. Applications in plant sciences, 8(8), e11383.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
## Not yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
