\name{BootstrapSmacof}
\alias{BootstrapSmacof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap on the distance matrices used for MDS with Smacof
}
\description{
Obtains bootstrap replicates of a distance matrix using ramdom samples or permuatations of a distance matrix. The object is to estimate the sampling variability of the results of the Smacof algorithm. 
}
\usage{
BootstrapSmacof(D, W=NULL, Model=c("Identity", "Ratio", "Interval", "Ordinal"), 
                dimsol=2, maxiter=100, maxerror=0.000001, StandardizeDisparities=TRUE,
                ShowIter=TRUE, nB=200, ProcrustesRot=TRUE, 
                method=c("Sampling", "Permutation"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
A distance matrix
}
  \item{W}{
A diagonal matrix containing waiths for the rows of D
}
  \item{Model}{
Mesurement level of the distances
}
  \item{dimsol}{
Dimension of the solution
}
  \item{maxiter}{
Maximum number of iterations for the smacof algorithm
}
  \item{maxerror}{
Tolerance for the smacof algorithm
}
  \item{StandardizeDisparities}{
Should the disparities be standardized in the smacof algorithm?
}
  \item{ShowIter}{
Should the information on each ieration be printed on the screen?
}
  \item{nB}{
Number of Bootstrap replications
}

  \item{ProcrustesRot}{
Should each replication be rotated to match the initial solution?
}
  \item{method}{
The replications are obtained "Sampling" or "Permutating" the residuals.
}

}
\details{
The function calculates bootstrap confidence intervals for coordinates and different stress measures using a distance matrix as a basis. The funcion uses random sampling or permutations of the residuals to obtain the bootstrap replications. The procedure preserves the length of the points in the multidimensional space perturbating only the angles among the vectors. It is done so to preserve the property of positiveness of the diagonal elements of the scalar product matrices. The procedure may result into a scalar product that does not have an euclidean configuration and then has some negative eigenvalues; to avoid this problem the negative eigenvalues are removed to approximate the perturbated matrix by the closest with the required properties.

It is well known that the eigenvectors of a matrix are unique except for reflections, that is, if we change the sign of each component of the eigenvector we have the same solution. If that happens, an unwanted increase in the variability due to this artifact may invalidate the results. To avoid this we can calculate the scalar product of each eigenvector of the initial matrix with the corresponding eigenvector of the bootstrap replicate and change the signs of the later if the result is negative.

Another artifact of the procedure may arise when the dimension of the solution is higher than 1 because the eigenvectors of a replicate may generate the same subspace although are not in the same directions, i. e., the subspace is referred to a different system. That also may produce an unwanted increase of the variability that invalidates the results. To avoid this, every replicate may be rotated to match as much as possible the subspace generated by the eigenvectors of the initial matrix. This is done by Procrustes Analysis, taking the rotated matrix as solution. 
The solution to this problem is also a sulution to the reflection, then only this problem is considered.
}
\value{
Returns an object of class "PCoABootstrap" with the information for each bootstrap replication.
  \item{Info}{Information about the procedure}
  \item{InitialDistance}{Initial distance}
  \item{RawStress}{A vector containing the raw stress for all the bootstrap replicates}
  \item{stress1}{A vector containing the value of the stress1 formula for all the bootstrap replicates}
  \item{stress2}{A vector containing the value of the stress2 formula for all the bootstrap replicates}
  \item{sstress1}{A vector containing the value of the sstress1 formula for all the bootstrap replicates}
  \item{sstress2}{A vector containing the value of the sstress2 formula for all the bootstrap replicates}
  \item{Coordinates}{A list with a component for each object. A component contains the coordinates of an object for all the bootstrap replicates (in columns)}
  \item{NReplicates}{Number of bootstrap replicates}
}
  
\references{
Efron, B.; Tibshirani, RJ. (1993). An introduction to the bootstrap. New York: Chapman and Hall. 436p.

Ringrose, T. J. (1992). Bootstrapping and correspondence analysis in archaeology. Journal of Archaeological Science, 19(6), 615-629.

MILAN, L., & WHITTAKER, J. (1995). Application of the parametric bootstrap to models that incorporate a singular value decomposition. Applied statistics, 44(1), 31-49.

Jacoby, W. G., & Armstrong, D. A. (2014). Bootstrap Confidence Regions for Multidimensional Scaling Solutions. American Journal of Political Science, 58(1), 264-278.
}
\author{
Jose L. Vicente-Villardon
}


\seealso{
\code{\link{BootstrapScalar}}
}
\examples{
\donttest{data(spiders)
D=BinaryProximities(spiders, coefficient="Jaccard", transformation="sqrt(1-S)")
DB=BootstrapDistance(D$Proximities)}
}

