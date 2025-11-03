test_that("liftOver", {

  region <- data.frame(seqnames = 1, start = 2000000, end = 2100000)

  expect_false((region %>% pull(start)) == (region %>% liftOverHg19() %>% as.data.frame() %>% pull(start)))
  expect_false((region %>% pull(start)) == (region %>% plyranges::as_granges() %>% liftOverHg19() %>% as.data.frame() %>% pull(start)))
})

test_that("Utility functions", {
  expect_no_error(exampleTumourNormal %>% makeTable() %>% qseaTableToChrGRanges())

  expect_equal(exampleTumourNormal %>% getWindowNames() %>% length(), exampleTumourNormal %>% getRegions() %>% length() )

  expect_equal(exampleTumourNormal %>% getPattern(), exampleTumourNormal@parameters$enrichmentPattern)

  expect_equal(exampleTumourNormal %>% convertToArrayBetaTable %>% dim(), c(724, 11) )

  expect_error(data.frame(a = 1, b = 2) %>% getSampleNames(), regexp = "only on a qseaSet")

})


test_that("calculateFractionReadsInGranges", {
  rgn <- data.frame(seqnames = 7, start = 25000000, end = 26000000)

  df0 <- exampleTumourNormal %>% calculateFractionReadsInGRanges(rgn, 0)
  df5 <- exampleTumourNormal %>% calculateFractionReadsInGRanges(rgn, 5)

  expect_equal(df0 %>% pull(sample_name), exampleTumourNormal %>% getSampleNames())
  expect_true(all(df5$afterOverBackNum <= df0$afterOverBackNum))
  expect_true(all(df5$initialOverBackNum <= df0$initialOverBackNum))
  expect_true(all(df5$initialOverBackNum != df0$initialOverBackNum))
})


test_that("getTables", {

  expect_equal(exampleTumourNormal %>% getCountTable() %>% nrow(), exampleTumourNormal %>% getRegions() %>% length())
  expect_equal(exampleTumourNormal %>% getNRPMTable() %>% nrow(), exampleTumourNormal %>% getRegions() %>% length())
  expect_equal(exampleTumourNormal %>% getDataTable() %>% nrow(), exampleTumourNormal %>% getRegions() %>% length())
  expect_equal(exampleTumourNormal %>% getBetaTable() %>% nrow(), exampleTumourNormal %>% getRegions() %>% length())

})

test_that("Downsampling", {

  randomSet <- qsea::getExampleQseaSet(repl = 2, expSamplingDepth = 10000)

  downsampleSet <- downSample(randomSet, nReads = 1000)
  expect_equal(downsampleSet %>% getCounts() %>% colSums() %>% unique(), 1000)
  expect_true(all(downsampleSet@count_matrix <= randomSet@count_matrix))

  downsampleSet <- downSample(randomSet, nReads = 100)
  expect_equal(downsampleSet %>% getCounts() %>% colSums() %>% unique(), 100)
  expect_true(all(downsampleSet@count_matrix <= randomSet@count_matrix))

  expect_error(downSample(randomSet, nReads = 100000))

})