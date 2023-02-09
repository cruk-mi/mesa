test_that("Annotation getting works", {

  expect_equal(getAnnotation(exampleTumourNormal, sampleAnnotation = tumour, useGroups = FALSE) %>% dim(),c(10,1))
  expect_equal(getAnnotation(exampleTumourNormal, sampleAnnotation = "tumour", useGroups = FALSE) %>% dim(), c(10,1))
  expect_equal(getAnnotation(exampleTumourNormal, sampleAnnotation = c("tumour","type"),  useGroups = FALSE)  %>% dim(), c(10,2))
  expect_equal(getAnnotation(exampleTumourNormal, sampleAnnotation = c(tumour,type), useGroups = FALSE)  %>% dim(), c(10,2))

  expect_equal(getAnnotation(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = tumour, useGroups = TRUE)  %>% dim(), c(2,1))
  expect_equal(getAnnotation(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = "tumour", useGroups = TRUE)  %>% dim(), c(2,1))

  expect_error(getAnnotation(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = c("tumour","type"), useGroups = TRUE)  %>% dim())
  expect_error(getAnnotation(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = c(tumour,type), useGroups = TRUE)  %>% dim())

  expect_no_error(plotCorrelationMatrix(exampleTumourNormal, sampleAnnotation = tumour, useGroups = FALSE))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal, sampleAnnotation = "tumour", useGroups = FALSE))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal, sampleAnnotation = c("tumour","type"),  useGroups = FALSE))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal, sampleAnnotation = c(tumour,type), useGroups = FALSE))

  expect_no_error(plotCorrelationMatrix(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = tumour, useGroups = TRUE))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = "tumour", useGroups = TRUE))

  expect_error(plotCorrelationMatrix(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = c("tumour","type"), useGroups = TRUE))
  expect_error(plotCorrelationMatrix(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = c(tumour,type), useGroups = TRUE))

  regions <- getRegions(exampleTumourNormal)[1:10]
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal, regionsToOverlap = regions, sampleAnnotation = tumour, useGroups = FALSE))
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal, regionsToOverlap = regions, sampleAnnotation = "tumour", useGroups = FALSE))
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal, regionsToOverlap = regions, sampleAnnotation = c("tumour","type"),  useGroups = FALSE))
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal, regionsToOverlap = regions, sampleAnnotation = c(tumour,type), useGroups = FALSE))

  expect_no_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = tumour, useGroups = TRUE))
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = "tumour", useGroups = TRUE))

  expect_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = c("tumour","type"), useGroups = TRUE))
  expect_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = c(tumour,type), useGroups = TRUE))
}
)

test_that("Testing hg38 related annotation/plotting functions", {

  expect_no_error(plotGeneHeatmap(exampleTumourNormal, gene = "HOXA10", sampleAnnotation = tumour, useGroups = FALSE))
  expect_no_error(plotGeneHeatmap(exampleTumourNormal, gene = "HOXA10", sampleAnnotation = "tumour", useGroups = FALSE))
  expect_no_error(plotGeneHeatmap(exampleTumourNormal, gene = "HOXA10", sampleAnnotation = c("tumour","type"),  useGroups = FALSE))
  expect_no_error(plotGeneHeatmap(exampleTumourNormal, gene = "HOXA10", sampleAnnotation = c(tumour,type), useGroups = FALSE))

  expect_no_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = tumour, useGroups = TRUE))
  expect_no_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = "tumour", useGroups = TRUE))

  expect_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = c("tumour","type"), useGroups = TRUE))
  expect_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = c(tumour,type), useGroups = TRUE))

  expect_no_error(exampleTumourNormal %>%
                    plotGeneHeatmap(gene = "HOXA10"))

  expect_no_error(exampleTumourNormal %>%
                    plotGeneHeatmap(gene = "HOXA10", normMethod = "nrpm"))

  expect_no_error(exampleTumourNormal %>%
                    plotGeneHeatmap(gene = "HOXA10", normMethod = "nrpm", maxScale = 3))

  expect_no_error(exampleTumourNormal %>%
                    plotGeneHeatmap(gene = "HOXA10", normMethod = "nrpm", maxScale = 3,
                                    sampleAnnotation = "tumour" ))

  expect_no_error(exampleTumourNormal %>%
                    plotCNVheatmap(tumour))

  expect_equal(exampleTumourNormal %>%
                 subsetWindowsBySignal(samples = "_N$",
                                       fn = max,
                                       threshold = 1,
                                       aboveThreshold = FALSE,
                                       normMethod = "nrpm") %>%
                 getRegions() %>%
                 length(), 601)
  
  expect_equal(exampleTumourNormal %>%
                 subsetWindowsBySignal(samples = c("Colon1_T","Colon2_T"),
                                       fn = max,
                                       threshold = 1,
                                       aboveThreshold = FALSE,
                                       normMethod = "nrpm") %>%
                 getRegions() %>%
                 length(), 604)
  
  expect_equal(exampleTumourNormal %>%
                 subsetWindowsBySignal(samples = c("Colon1_T","Colon2_T"),
                                       fn = max,
                                       threshold = 1,
                                       aboveThreshold = TRUE,
                                       normMethod = "nrpm") %>%
                 getRegions() %>%
                 length(), 215)

  expect_equal(exampleTumourNormal %>%
                 summariseAcrossWindows(windowsToUse = getRegions(.)) %>%
                 nrow(), exampleTumourNormal %>% getSampleNames() %>% length())

  expect_equal(exampleTumourNormal %>% pull(group) %>% length(), exampleTumourNormal %>% getSampleNames() %>% length())
  expect_equal(exampleTumourNormal %>% pull(tumour) %>% unique(), c("Normal", "Tumour"))
  expect_equal(exampleTumourNormal %>% pull(type) %>% unique() %>% length(), 5)

})

test_that("Testing general functionality", {
  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000000) %>%
    mutate(patient = stringr::str_remove(sample_name,"[TN]$"),
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

  expect_no_error(DMRs <- exampleTumourNormal %>%
    calculateDMRs(variable = "type",
                  contrasts = tibble::tibble(sample1 = c("LUAD","LUSC","CRC"),
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
                 summariseDMRsByGene() %>%
                 dim(),c(28,8))

})

test_that("Multiple DMRs", {
  expect_no_error(DMRs <- exampleTumourNormal %>%
    calculateDMRs(variable = "type",
                  contrasts = "all"))
})
