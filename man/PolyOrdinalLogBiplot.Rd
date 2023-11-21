\name{PolyOrdinalLogBiplot}
\alias{PolyOrdinalLogBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Factor Analysis Biplot based on polychoric correlations
}
\description{
Calculates a biplot for ordinal data based on polychoric correlations
}
\usage{
PolyOrdinalLogBiplot(X, dimension = 3, method = "principal", 
rotate = "varimax", RescaleCoordinates = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A matrix of ordinal data
}
  \item{dimension}{
Number of dimensiona to retain
}
  \item{method}{
Principal components (principal) or factor analysis (fa)
}
  \item{rotate}{
Rotation for the analysis
}
  \item{RescaleCoordinates}{
Rescale coordinates as in a continuous data biplot
}
  \item{\dots}{
Any aditional arguments for the principal and fa functions
}
}
\details{
The procedure calculates 
}
\value{
A biplot (Continuous or ordinal)
}
\author{
Jose Luis Vicente Villardon
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fa}}, \code{\link{principal}}
}
\examples{
## Not Yet
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
