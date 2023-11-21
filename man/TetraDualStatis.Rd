\name{TetraDualStatis}
\alias{TetraDualStatis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dual STATIS-ACT for binary data based on Tetrachoric Correlations
}
\description{
Dual STATIS-ACT for binary data based on Tetrachoric Correlations
}
\usage{
TetraDualStatis(X, dimens = 2, SameInd = FALSE, RotVarimax = FALSE, 
               OptimMethod = "L-BFGS-B", penalization = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A three way binary data matrix
}
  \item{dimens}{
Dimension of the solution
}
  \item{SameInd}{
Are the individuals the same in all occassions?
}
  \item{RotVarimax}{
Should the solution be rotated?
}
  \item{OptimMethod}{
Optimization method for the gradients
}
  \item{penalization}{
Penalization for the ridge solution
}
}
\details{
The general aim of STATIS-ACT methods is to extract information common to a set of datasets with the same individuals. They will also be represented as a Euclidean configuration or map of points (or vectors), in the same way as in Principal Component Analysis (PCA) or Principal Coordinate Analysis (PCoA). If the object is to analyze the variables and the correlation structures between them we will use a Factor Analysis (FA). When we have tables in which we measure a set of common variables and we want to obtain a consensus structure of all of them, we will use the named STATIS-Dual.

The method was initially designed to work with individuals common to all the tables, but in this work, we will focus on the dual version, which works with variables common to all of them.

When we have several tables of binary dataset, the classical methods for continuous data are not suitable. If the individuals are the same in all tables, we can use a STATIS based on distances, also known as DISTATIS. El procedimiento consiste en calcular una matriz de distancias a partir de para un coeficiente de similaridad para datos binarios. Las distancias se convierten en productos escalares, como en ACoP, y se trabaja a partir de ellos como en el STATIS tradicional.

When we have common variables, and we are interested in the association between them, we could use a coefficient that, instead of similarity, shows the association between the variables. In this work we propose the use of the tetrachoric correlation matrix for each table and develop the necessary adaptations to the method.

}
\value{
An object with the results
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Laura Vicente-Gonzalez, José Luis Vicente-Villardon
}

\examples{
# Not yet
}
