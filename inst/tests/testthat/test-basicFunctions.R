test_that("plotColors works", {
  colors <- plotColors(x = RepSeqData)
  expect_equal(any(is.na(colors)), FALSE)
})

test_that("countFeatures works", {
level_statistics <- countFeatures(x = RepSeqData,
                                  level = "J",
                                  group=c("sex", "F"),
                                  scale="count")
expect_equal(dim(level_statistics), c(50,5))
})