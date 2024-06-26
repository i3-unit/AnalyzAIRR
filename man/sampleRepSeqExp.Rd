% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diversityFunctions.R
\name{sampleRepSeqExp}
\alias{sampleRepSeqExp}
\title{Down-sampling of repertoires}
\usage{
sampleRepSeqExp(
  x,
  sample.size = min(mData(x)$nSequences),
  rngseed = FALSE,
  replace = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{x}{an object of class \code{\linkS4class{RepSeqExperiment}}}

\item{sample.size}{an integer indicating the desired down-sampled size.
The default is the smallest repertoire size among all samples of the dataset.}

\item{rngseed}{a integer used as seed for a reproducible result.}

\item{replace}{a boolean indicating if the resampling should be performed
with or without replacement. Default is TRUE.}

\item{verbose}{a boolean indicating whether or not to show the details of
every computation step within the function. Default is TRUE.}
}
\value{
a new \code{\linkS4class{RepSeqExperiment}}  with the downsized data.
}
\description{
This function downsamples all repertoires in the dataset.

Users can choose the value to which all the samples are down-sampled.
If not specified, the lowest number of sequences across all samples within
the dataset will be used.

This strategy can be applied when studying different cell subsets with
significant differences in their repertoire sizes.
}
\examples{

data(RepSeqData)

RepSeqData_ds<- sampleRepSeqExp(x = RepSeqData,
                                rngseed = 1234,
                                replace = TRUE)

RepSeqData_ds<- sampleRepSeqExp(x = RepSeqData,
                                rngseed = FALSE,
                                replace = FALSE)

}
