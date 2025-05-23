% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadData.R
\name{RepSeqExp}
\alias{RepSeqExp}
\title{An alternative method for the building of a RepSeqExperiment object}
\usage{
RepSeqExp(clonotypetab, sampleinfo = NULL)
}
\arguments{
\item{clonotypetab}{a single clonotype table containing the previously cited columns}

\item{sampleinfo}{a data frame containing:

- a column with the sample names. Names should match the base names of the corresponding files and their order. This column should be assigned as row.names when the metadata file is loaded. See the example below.

- any additional columns with relevant information for the subsequent analyses. Group columns must be transformed into factors after loading. See the example below.

No specific column names are required.}
}
\value{
an object of class \code{RepSeqExperiment} that is  used in all the analytical metrics proposed by the AnalyzAIRR package. See \code{\link{RepSeqExperiment-class}} for more details.
}
\description{
This function can be used to build a RepSeqExperiment object using aligned files that were not produced by a AnalyzAIRR-supported aligning tool. The loaded clonotype table must however contain the following required columns:

 - sample_id: sample names

 - V: Variable gene name

 - J: Joining gene name

 - aaCDR3: amino acid CDR3 sequence

 - ntCDR3: nucleotide CDR3 sequence

 - aaClone: Full clone sequence including the V gene, the amino acid CDR3 sequence and the J gene
 
 - VJ: V-J gene combination using V and J gene names

 - ntClone: Full clone sequence including the V gene, the nucleotide CDR3 sequence and the J gene

 - count: the occurrence of the clone

Clonotype tables must only contain a single chain. No paired-chain analysis are provided by the AnalyzAIRR package.
Pre-filtered files obtained using the \code{\link{filterClones}} function can be used as input.
}
\keyword{internal}
