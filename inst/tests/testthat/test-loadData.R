test_that("multiplication works", {
  l <-   list.files(system.file(file.path('extdata/mixcr'),
                                 package = 'AnalyzAIRR'),
                                 full.names = TRUE)

  
  metaData <- read.table(system.file(file.path('extdata/sampledata.txt'),
                                     package = 'AnalyzAIRR'),sep = "\t",row.names = 1, header = TRUE)
  metaData$cell_subset <- factor(metaData$cell_subset)
  metaData$sex <- factor(metaData$sex)
  
  dataset <- AnalyzAIRR::readAIRRSet(fileList = l,
                         fileFormat = "MiXCR",
                         chain = "TRA",
                         sampleinfo = metaData,
                         filter.singletons = FALSE,
                         aa.th=9,
                         outFiltered = FALSE,
                         raretab = FALSE,
                         cores=1L)
  expect_equal(is.RepSeqExperiment(dataset), TRUE)
})
