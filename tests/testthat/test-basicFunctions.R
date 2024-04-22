test_that("plotColors works", {
  colors <- AnalyzAIRR::plotColors(x = RepSeqData)
  expect_equal(any(is.na(colors)), FALSE)
})

test_that("countFeatures works", {
level_statistics <- AnalyzAIRR::countFeatures(x = RepSeqData,
                                  level = "J",
                                  group=c("sex", "F"),
                                  scale="count")
expect_equal(dim(level_statistics), c(44,5))
})