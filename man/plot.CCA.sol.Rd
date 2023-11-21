\name{plot.CCA.sol}
\alias{plot.CCA.sol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the solution of a Canonical Correspondence Analysisis
}
\description{
Plots the solution of a Canonical Correspondence Analysisis using similar parameters to the continuous biplot
}
\usage{
\method{plot}{CCA.sol}(x, A1 = 1, A2 = 2, ShowAxis = FALSE, margin = 0,
                 PlotSites = TRUE, PlotSpecies = TRUE, PlotEnv = TRUE,
                 LabelSites = TRUE, LabelSpecies = TRUE, LabelEnv =
                 TRUE, TypeSites = "wa", SpeciesQuality = FALSE,
                 MinQualityVars = 0.3, dp = 0, pr = 0, PlotAxis =
                 FALSE, TypeScale = "Complete", ValuesScale =
                 "Original", mode = "a", CexSites = NULL, CexSpecies =
                 NULL, CexVar = NULL, ColorSites = NULL, ColorSpecies =
                 NULL, ColorVar = NULL, PchSites = NULL, PchSpecies =
                 NULL, PchVar = NULL, SizeQualSites = FALSE,
                 SizeQualSpecies = FALSE, SizeQualVars = FALSE,
                 ColorQualSites = FALSE, ColorQualSpecies = FALSE,
                 ColorQualVars = FALSE, SmartLabels = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The results of a CCA model
}
  \item{A1}{
Dimension for the first axis
}
  \item{A2}{
Dimension for the second axis
}
  \item{ShowAxis}{
Logical variable to control if the coordinate axes should appear in the plot. The default value is FALSE because for most of the biplots its presence is irrelevant.
}
  \item{margin}{
Margin for the labels in some of the biplot modes (percentage of the plot width). Default is 0. Increase the value if the labels are not completely plotted.
}
  \item{PlotSites}{
Should the sites be plotted?
}
  \item{PlotSpecies}{
Should the species be plotted?
}
  \item{PlotEnv}{
Should the environmental variables be plotted?
}
  \item{LabelSites}{
Labels for the sites
}
  \item{LabelSpecies}{
Labels for the species
}
  \item{LabelEnv}{
Labels for the environmental variables.
}
  \item{TypeSites}{
Type for the sites plot
}
  \item{SpeciesQuality}{
Quality for the species
}
  \item{MinQualityVars}{
Minimum quality to plot a variable
}
  \item{dp}{
A set of indices with the variables that will show the projections of the individuals.
}
  \item{pr}{
A set of indices with the individuals to show the projections on the variables.
}
  \item{PlotAxis}{
Should the axis be plotted?
}
  \item{TypeScale}{
Type of scale to use : "Complete", "StdDev" or "BoxPlot"
}
  \item{ValuesScale}{
Values to show on the scale: "Original" or "Transformed"
}
  \item{mode}{
Mode of the biplot: "p", "a", "b", "h", "ah" and "s".
}
  \item{CexSites}{
Size for the symbols and labels of the sites. Can be a single common size for all the points or a vector with individual sizes.
}
  \item{CexSpecies}{
Size for the symbols and labels of the species. Can be a single common size for all the points or a vector with individual sizes.
}
  \item{CexVar}{
Size for the symbols and labels of the variables. Can be a single common size for all the points or a vector with individual sizes.
}
  \item{ColorSites}{
Color for the symbols and labels of the sites. Can be a single common color for all the points or a vector with individual colors.
}
  \item{ColorSpecies}{
Color for the symbols and labels of the species. Can be a single common color for all the points or a vector with individual colors.
}
  \item{ColorVar}{
Color for the symbols and labels of the variables. Can be a single common color for all the points or a vector with individual colors.
}
  \item{PchSites}{
Symbol for the sites points. See \code{help(points)} for details.
}
  \item{PchSpecies}{
Symbol for the species points. See \code{help(points)} for details.
}
  \item{PchVar}{
Symbol for the variables points. See \code{help(points)} for details.
}
  \item{SizeQualSites}{
Should the size of the site points be related to their qualities of representation (predictiveness)?
}
  \item{SizeQualSpecies}{
Should the size of the species points be related to their qualities of representation (predictiveness)?
}
  \item{SizeQualVars}{
Should the size of the variables points be related to their qualities of representation (predictiveness)?
}
  \item{ColorQualSites}{
Should the color of the sites points be related to their qualities of representation (predictiveness)?
}
  \item{ColorQualSpecies}{
Should the color of the species points be related to their qualities of representation (predictiveness)?
}
  \item{ColorQualVars}{
Should the color of the variables points be related to their qualities of representation (predictiveness)?
}
  \item{SmartLabels}{
Plot the labels in a smart way
}
  \item{...}{
Aditional graphical parameters.
}
}
\details{
The plotting procedure is similar to the one used for continuous biplots including the calibration of the environmental variables.
}
\value{
No value returned
}
\references{
CCA
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
 \code{\link{plot.ContinuousBiplot}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
}

