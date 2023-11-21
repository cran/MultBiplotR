\name{BinaryPLSFit}
\alias{BinaryPLSFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary PLS Regression.
}
\description{
Fits Binary PLS regression.
}
\usage{
BinaryPLSFit(Y, X, S = 2, tolerance = 5e-06, maxiter = 100, show = FALSE, 
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
Optimization methods from optimr
}
  \item{seed}{
Seed. By default is 0.
}
}
\details{
Fits Binary PLS Regression. It is used for a higher level function.
}
\value{
The PLS fit used by the BinaryPLSR function.
}
\references{
Ugarte Fajardo, J., Bayona Andrade, O., Criollo Bonilla, R., Cevallos‐Cevallos, J., Mariduena‐Zavala, M., Ochoa Donoso, D., & Vicente Villardon, J. L. (2020). Early detection of black Sigatoka in banana leaves using hyperspectral images. Applications in plant sciences, 8(8), e11383.

Vicente-Gonzalez, L., & Vicente-Villardon, J. L. (2022). Partial Least Squares Regression for Binary Responses and Its Associated Biplot Representation. Mathematics, 10(15), 2580.
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
