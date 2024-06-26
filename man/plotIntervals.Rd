% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{plotIntervals}
\alias{plotIntervals}
\title{Visualization of the clonal distribution per interval}
\usage{
plotIntervals(
  x,
  level = c("aaClone", "ntClone", "ntCDR3", "aaCDR3"),
  grouped = FALSE,
  colorBy = NULL,
  facetBy = NULL,
  label_colors = NULL,
  show_stats = FALSE,
  fractions = c("count", "frequency")
)
}
\arguments{
\item{x}{an object of class  \code{\linkS4class{RepSeqExperiment}}}

\item{level}{a character specifying the level of the repertoire to be taken into account when calculating the clonal distribution. Should be one of "aaClone","ntClone", "ntCDR3" or "aaCDR3".}

\item{grouped}{a character indicating one or multiple groups to be compared. A Wilcoxon test is thus performed and adjusted p-values using the Holm method are shown. Colors are attributed to the different groups within the first column, and a facet is applied on the second column. If not specified, no statistical tests will be performed, and calculated values for each sample_id will be represented.}

\item{colorBy}{a character indicating a column name in mData. Colors are thus attributed to the different groups within this column. The chosen column must be of class factor.}

\item{facetBy}{a vector of character indicating one or two column names in mData to apply a facet on.}

\item{label_colors}{a list of colors for each variable in ColorBy. See \code{\link{plotColors}}. If NULL, default colors are used.}

\item{show_stats}{whether to statistically compare groups}

\item{fractions}{whether intervals should be determined in count or frequency}
}
\description{
This function plots two histograms of the clonal distribution per a set of intervals in all the samples within the dataset. 

The plot titled "Distribution" calculates the proportion of each interval in the whole repertoire, whereas the one titled "Cumulative frequency" shows the cumulative frequency of the sequences within each interval.

This could allow a global view of the repertoire fraction contributing the most to the repertoire. For instance, top sequences belonging to the highest interval often constitute a low fraction of the whole repertoire but contribute more significantly in terms of cumulative frequency in view of their high occurrence.

Samples can be statistically compared in each interval using the \code{groupBy} parameter.
}
\examples{

data(RepSeqData)

plotIntervals(x = RepSeqData, level="aaCDR3", facetBy="cell_subset", fractions="count")

plotIntervals(x = RepSeqData, level="ntCDR3", colorBy="cell_subset", grouped=TRUE, 
              fractions="frequency", show_stats=TRUE )

}
