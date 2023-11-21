\name{t3pcovr}
\alias{t3pcovr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tucker 3 Principal Covariates Regression
}
\description{
Tucker 3 Principal Covariates Regression
}
\usage{
t3pcovr(X, Y, I, J, K, L, r1 = 2, r2 = 2, r3 = 2, 
       conv = 1e-06, OriginalAlfa = 0.5, AlternativeLossF = 1, 
       nRuns = 100, StartSeed = 0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{X}{A two way data matrix with the predictors.}
  \item{Y}{A three way data matrix with the responses.}
  \item{I}{Number of elements of first mode of 3D/2D (the common mode: rows)}
  \item{J}{ number of elements of second mode of 3D (columns 3D) }
   \item{K}{number of elements of third mode of 3D (slabs)}
  \item{L}{number of elements of second mode of 2D (columns 2D)}
  \item{r1}{Number of extracted components for the A-mode}
  \item{r2}{Number of extracted components for the B-mode}
  \item{r3}{Number of extracted components for the C-mode}
  \item{conv}{value for convergence (tolerance value)}
  \item{OriginalAlfa}{(0-1): importance that degree reduction and prediction have in the analysis}
  \item{AlternativeLossF}{Using the alternative loss function? 0 = no (use original loss function: weighted SSQ;               weighted met alfa) 1 = yes (use weighted loss function with scaled SSQ: scaled by the SSQ in X and y )}
  \item{nRuns}{Number of runs}
  \item{StartSeed}{Seed for the analysis}
}

\details{
In behavioral research it is very common to have to deal with  several data sets which include information relative to the same set of individuals, in such a way that one data set tries to explain the others. The class of models known as PCovR focuses on the analysis of a three-way data array explained by a two-way data matrix. In this paper the Tucker3-PCovR model is proposed that is a particular case of PCovR class. Tucker3-PCovR model reduces the predictors to a few components and predict the criterion by using these components and, at the same time, the three way data is fitted through the Tucker3 model. Both, the reduction of the predictors and the prediction of the criterion are done simultaneously. An alternating least squares algorithm to estimate the Tucker3-PCovR model is proposed. A biplot representation to facilitate the interpretation of the results is presented. A couple of applications are made to coupled empirical data sets related to the field of psychology. 
}
\value{
  \item{A}{Component matrix for the A-mode)}
  \item{B1}{Component matrix for the B-mode}
  \item{C}{Component matrix for the C-mode}
  \item{H}{Matrized core array (frontal slices)}
  \item{B2}{Loading matrix of components (components x predictors)}
  \item{...}{Further arguments}
}


\references{\cite{De Jong, S., & Kiers, H. A. (1992). Principal covariates regression: Part I. Theory. Chemometrics and Intelligent Laboratory Systems , 155-164.}

\cite{Marlies Vervloet, Henk A. Kiers, Wim Van den Noortgate, Eva Ceulemans (2015). PCovR: An R Package for Principal Covariates Regression. Journal of Statistical Software, 65(8), 1-14. URL http://www.jstatsoft.org/v65/i08/.}

\cite{Smilde, A. K., Bro, R., & Geladi, P. (2004). Multi-way analysis with applications in the chemical sciences. Chichester, UK: Wiley.}}

\author{
Elisa Frutos Bernal (\email{efb@usal.es})
}

\examples{
#Not yet
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
