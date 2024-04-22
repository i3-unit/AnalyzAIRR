test_that("RepSeqData works", {
  res<- AnalyzAIRR::diversityIndices(RepSeqData, level="V")
  expect_equal(sum(res$shannon),32.27)
})
