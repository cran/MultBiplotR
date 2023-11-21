\name{plot.Unfolding}
\alias{plot.Unfolding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots an Unfolding Representation
}
\description{
Plots an Unfolding Representation
}
\usage{
\method{plot}{Unfolding}(x, A1 = 1, A2 = 2, ShowAxis = FALSE,
margin = 0.1, PlotSites = TRUE, PlotSpecies = TRUE, PlotEnv = TRUE,
LabelSites = TRUE, LabelSpecies = TRUE, LabelEnv = TRUE, 
SpeciesQuality = FALSE, MinQualityVars = 0, dp = 0, 
PlotAxis = FALSE, TypeScale = "Complete", ValuesScale = "Original", 
mode = "h", CexSites = NULL, CexSpecies = NULL, CexVar = NULL, 
ColorSites = NULL, ColorSpecies = NULL, ColorVar = NULL, 
PchSites = NULL, PchSpecies = NULL, PchVar = NULL, 
SizeQualSites = FALSE, SizeQualSpecies = FALSE, 
SizeQualVars = FALSE, ColorQualSites = FALSE, 
ColorQualSpecies = FALSE, ColorQualVars = FALSE, SmartLabels = FALSE, 
PlotTol = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class Unfolding
}
  \item{A1}{
Axis 1 of the representation.
}
  \item{A2}{
Axis 1 of the representation.
}
  \item{ShowAxis}{
Should the axis be shown?
}
  \item{margin}{
Margin for the plot (precentage)
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
Should the sites be labelled?
}
  \item{LabelSpecies}{
Should the species be labelled?
}
  \item{LabelEnv}{
Should the environmental variables be labelled?
}
  \item{SpeciesQuality}{
Min species quality to plot
}
  \item{MinQualityVars}{
Minimum quality of a var to be plotted.
}
  \item{dp}{
A set of indices with the variables that will show the projections of the individuals.
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
  \item{PlotTol}{
Should the tolerances be plotted
}
  \item{...}{
Aditional graphical parameters.
}
}
\details{
Plots an Unfolding Representation
}
\value{
A plot of the unfolding representation.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
de Leeuw, J. (2005). Multidimensional unfolding. Encyclopedia of statistics in behavioral science.
}
\author{
Jose Luis Vicente-Villardon
}

\examples{
# Not yet
}
