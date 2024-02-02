test_that("multiplication works", {
  l <- list.files(file.path('../../inst/extdata/mixcr'),
                  full.names = TRUE)
  
  metaData <- read.table(file.path('../../inst/extdata/sampledata.txt'),sep = "\t",row.names = 1, header = TRUE)
  metaData$cell_subset <- factor(metaData$cell_subset)
  metaData$sex <- factor(metaData$sex)
  
  dataset <- readAIRRSet(fileList = l,
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
