\name{BiplotFPCA}
\alias{BiplotFPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
External Biplot for functional data from a functional PCA object.
}
\description{
The function calculates a biplot from a functional PCA object and the data used tocalculate it.
}
\usage{
BiplotFPCA(FPCA, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FPCA}{
Functional PCA object
}
  \item{X}{
Data used to calculate the fuctional PCA
}
}
\details{
The function calculates a biplot from a functional PCA object and the data used tocalculate it. At this moment the function calculates only an external biplot by regressing X o the funcional components. Furure versions will include the internal biplot.
}
\value{
A Continuous biplot object
}

\author{
José Luis Vicente Villardón
}


\examples{
# not yet
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
