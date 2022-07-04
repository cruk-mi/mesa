test_that("Testing hg38 related annotation/plotting functions", {
  expect_no_error(examplePairedTumourQset %>%
                    plotGeneHeatmap(gene = "HOXA10"))

  expect_no_error(examplePairedTumourQset %>%
                    plotGeneHeatmap(gene = "HOXA10", normMethod = "nrpm"))

  expect_no_error(examplePairedTumourQset %>%
                    plotGeneHeatmap(gene = "HOXA10", normMethod = "nrpm", maxScale = 3))

  expect_no_error(examplePairedTumourQset %>%
                    plotGeneHeatmap(gene = "HOXA10", normMethod = "nrpm", maxScale = 3,
                                    annotationCol = getAnnotationDataFrame(., "tumour") ))

  expect_no_error(examplePairedTumourQset %>%
                    plotCNVheatmap(tumour))

  expect_equal(examplePairedTumourQset %>%
                 removeExpressedWindows(samplesToFilterOut = "_N$",
                                        maxValue = 1,
                                        normMethod = "nrpm") %>%
                 getRegions() %>%
                 length(), 369)

  expect_equal(examplePairedTumourQset %>%
                 removeExpressedWindows(samplesToFilterOut = c("Colon1_T","Colon2_T"),
                                        maxValue = 1,
                                        normMethod = "nrpm") %>%
                 getRegions() %>%
                 length(), 264)

  expect_equal(examplePairedTumourQset %>%
                 keepExpressedWindows(samplesToFilterOn = c("Colon1_T","Colon2_T"),
                                        minValue = 1,
                                        normMethod = "nrpm") %>%
                 getRegions() %>%
                 length(), 555)

  expect_equal(examplePairedTumourQset %>%
                 getNormalisedReadSum(GRanges = getRegions(.)) %>%
                 nrow(), examplePairedTumourQset %>% getSampleNames() %>% length())

})

test_that("Testing general functionality", {
  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000000) %>%
    mutateQset(patient = stringr::str_remove(sample_name,"[TN]$"),
               variableWithOneLevel = "Test",
               experiment = ifelse( stringr::str_detect(sample_name,"[1234]"),"A","B"),
               experimentConfounded = ifelse( stringr::str_detect(sample_name,"[1234]N"),"A","B"),
    )

  expect_equal(randomSet %>% getBetaTable() %>% nrow(), randomSet %>% getRegions() %>% length())
  expect_equal(randomSet %>% getBetaTable() %>% dplyr::select(matches("Sim")) %>% ncol(), randomSet %>% getSampleNames() %>% length())
  expect_equal(randomSet %>% getBetaTable(groupMeans = TRUE) %>% dplyr::select(matches("Tum|Norm")) %>% ncol(), randomSet %>% getSampleGroups() %>% length())

  expect_equal(randomSet %>% getNRPMTable() %>% nrow(), randomSet %>% getRegions() %>% length())
  expect_equal(randomSet %>% getNRPMTable() %>% dplyr::select(matches("Sim")) %>% ncol(), randomSet %>% getSampleNames() %>% length())
  expect_equal(randomSet %>% getNRPMTable(groupMeans = TRUE) %>% dplyr::select(matches("Tum|Norm")) %>% ncol(), randomSet %>% getSampleGroups() %>% length())


})


test_that("Analysing DMRs", {

  expect_no_error(DMRs <- examplePairedTumourQset %>%
    calculateDMRs(variable = "type",
                  contrastsToDo = tibble::tibble(sample1 = c("LUAD","LUSC","CRC"),
                                                 sample2 = c("NormalLung","NormalLung","NormalColon"))))

  expect_equal(DMRs %>%
                 pivotDMRsLonger() %>%
                 nrow(), 81)

  expect_equal(DMRs %>%
                 summariseDMRsByContrast() %>%
                 dim(),c(3,4))

  expect_equal(DMRs %>%
                 pivotDMRsLonger() %>%
                 summariseDMRsByContrast() %>%
                 dim(),c(3,4))

  expect_equal(DMRs %>%
                 pivotDMRsLonger() %>%
                 summariseByGene() %>%
                 dim(),c(28,8))

})

test_that("Multiple DMRs", {
  DMRs <- examplePairedTumourQset %>%
    calculateDMRs(variable = "type",
                  contrastsToDo = "all")
})
