testthat::test_that("RepSeqData works", {
  res<- AnalyzAIRR::diversityIndices(RepSeqData, level="V")
  expect_equal(res$shannon[1],4.042)
})
