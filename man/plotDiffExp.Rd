% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{plotDiffExp}
\alias{plotDiffExp}
\title{Visualization of differential expression in a volcano plot}
\usage{
plotDiffExp(
  x,
  level = c("aaClone", "ntClone", "V", "J", "VJ", "ntCDR3", "aaCDR3"),
  group = c("cell_subset", "amTreg", "nTreg"),
  FC.TH = 2,
  PV.TH = 0.05,
  top = 10
)
}
\arguments{
\item{x}{an object of class  \code{\linkS4class{RepSeqExperiment}}}

\item{level}{a character specifying the level of the repertoire on which the diversity should be estimated. Should be one of "aaClone","ntClone", "V", "J", "VJ", "ntCDR3" or "aaCDR3".}

\item{group}{a vector of character indicating the column name in the mData slot, as well as the two groups to be compared.}

\item{FC.TH}{an integer indicating the log2FoldChange threshold. Default is 2.}

\item{PV.TH}{an integer indicating the adjusted pvalue threshold. Default is 0.05.}

\item{top}{an integer indicating the top n significant labels to be shown on the volcano plot. Default is 10.}
}
\description{
This function plots differentially expressed repertoire levels calculated using the \code{\link{diffExpGroup}} in a volcano plot.
}
\examples{

plotDiffExp(x = RepSeqData,
            level = "V",
            group = c("cell_subset", "amTreg", "nTreg"),
            top = 10,
            FC.TH = 1,
            PV.TH = 0.05)
            
}
