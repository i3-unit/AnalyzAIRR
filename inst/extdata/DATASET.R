## code to prepare `RepSeqData` dataset goes here
library(AnalyzAIRR)
# load sample information
setwd("/mnt/mukkuri/RepSeq/RS_Analysis/GPI/AnalyzAIRR/package")
sampleData <- read.table("/mnt/mukkuri/RepSeq/RS_Analysis/GPI/AnalyzAIRR/extdata/sampledata.txt", sep = "\t", header = TRUE, row.names=1)
sampleData$sex<- factor(sampleData$sex)
sampleData$cell_subset<- factor(sampleData$cell_subset)
# load aligner outputs
l <- list.files("/mnt/mukkuri/RepSeq/RS_Analysis/GPI/AnalyzAIRR/extdata/mixcr",pattern = ".txt", full.names=TRUE)
# create RepSeqExperiment
RepSeqData <- AnalyzAIRR::readAIRRSet(l,
                                 cores = 2L,
                                 fileFormat = "MiXCR",
                                 chain = "TRA",
                                 sampleinfo = sampleData,
                                 keep.ambiguous = FALSE,
                                 keep.unproductive = FALSE,
                                 filter.singletons=TRUE,
                                 aa.th = 8,
                                 outFiltered = TRUE,
                                 raretab = FALSE)
# use
usethis::use_data(RepSeqData,  compress="bzip2", overwrite = TRUE)
