% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadData.R
\name{parseMiXCR}
\alias{parseMiXCR}
\title{parse MiXCR output}
\usage{
parseMiXCR(path, chain = c("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"))
}
\arguments{
\item{path}{full path to the aligned file. Files can be loaded as gzipped.}

\item{chain}{character, the TCR chain to be analyzed. One of \code{A} or \code{B}. Default value is \code{A}.}
}
\value{
a data.table of 9 columns: \code{sample_id} name, \code{V} V gene, \code{J} J gene,\code{VJ} V-J gene combinations, \code{aaCDR3} amino acid CDR3 sequence, \code{ntCDR3} nucleotide CDR3 sequence, \code{ntClone} full clone sequence in nucleotide, \code{aaClone} full clone sequence in amino acid, \code{count} clone count. Clones are eliminated if a STOP codon (*) is detected in aaCDR3 chain, if the ntCDR3 length is not divisible by 3 or if the ntCDR3 sequence is ambiguous (contains a "N" base).
}
\description{
parse clonotype tables from MiXCR
}
\details{
function imports clonotype tables produced by the MiXCR aligner.
}
\keyword{internal}
