\name{plot.TetraDualStatis}
\alias{plot.TetraDualStatis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots an object of class "tetraDualStatis".
}
\description{
Plots an object the results of TetraDualStatis.
}
\usage{
\method{plot}{TetraDualStatis}(x, A1 = 1, A2 = 2, PlotType = "InterStructure", 
                    PlotRowTraj = FALSE, PlotVarTraj = FALSE, LabelTraj = "Begining", 
                    VarColorType = "Biplot", VarColors = NULL, VarLabels = NULL, 
                    RowColors = NULL, TableColors = NULL, RowRandomColors = FALSE, 
                    TypeTraj = "line", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class TetraDualStatis
}
  \item{A1}{
Dimension for the X-axis
}
  \item{A2}{
Dimension for the Y-axis
}
  \item{PlotType}{
Type of plot: "Biplot", "Compromise", "Correlations", "Contributions", "InterStructure".
}
  \item{PlotRowTraj}{
Should the row trajectories be plotted?
}
  \item{PlotVarTraj}{
Should the variables trajectories be plotted?
}
  \item{LabelTraj}{
Should the trajectories be labelled.
}
  \item{VarColorType}{
One of the following: "Biplot", "ByTable", "ByVar".
}
  \item{VarColors}{
User colors for the variables.
}
  \item{VarLabels}{
User labels for the variables.
}
  \item{RowColors}{
User colors for the rows.
}
  \item{TableColors}{
User colors for the different tables.
}
  \item{RowRandomColors}{
Should use random colors for the rows?
}
  \item{TypeTraj}{
Type of trajectory. Normally a line.
}
  \item{\dots}{
Additional graphical arguments.
}
}
\details{
Plots an object the results of TetraDualStatis.
}
\value{
The plot of the results
}

\author{
Laura Vicente-Gonzalez, Jose Luis Vicente-Villardon
}

\examples{
##---- Should be DIRECTLY executable !! ----

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
