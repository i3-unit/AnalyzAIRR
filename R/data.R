#' @title RepSeqData
#' @description A RepSeqExperiment object built using a published TCR repertoire dataset of healthy murine cell populations (Mhanna et al., 2021). Files aligned using MiXCR following the MiAIRR standardsare available in extdata/MiAIRR.
#' @format RepSeqExperiment object and gz-compressed txt/csv tab-delimited:
#' \describe{
#'   \item{\code{system.file(file.path('extdata/MiAIRR'), package='AnalyzAIRR')}}{Aligned data sets using MiXCR and exported under the MiAIRR standards.}
#'   \item{\code{system.file(file.path('extdata/sampledata.txt'), package='AnalyzAIRR')}}{Metadata}
#'}
#' @usage data(RepSeqData)
"RepSeqData"
