testthat::test_that("multiplication works", {
  l <-   list.files(system.file(file.path('extdata/MiAIRR'),
                                 package = 'AnalyzAIRR'),
                                 full.names = TRUE)

  
  metaData <- read.table(system.file(file.path('extdata/sampledata.txt'),
                                     package = 'AnalyzAIRR'),sep = "\t",row.names = 1, header = TRUE)
  metaData$cell_subset <- factor(metaData$cell_subset)
  metaData$sex <- factor(metaData$sex)
  
  dataset <- AnalyzAIRR::readAIRRSet(fileList = l,
                         fileFormat = "MiAIRR",
                         chain = "TRA",
                         sampleinfo = metaData,
                         filter.singletons = TRUE,
                         aa.th=9,
                         outFiltered = FALSE,
                         cores=1L)
  expect_equal(is.RepSeqExperiment(dataset), TRUE)
})
