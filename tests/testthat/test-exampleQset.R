test_that("Annotation getting works", {

  expect_equal(getAnnotation(exampleTumourNormal, sampleAnnotation = tumour, useGroupMeans = FALSE) %>% dim(),c(10,1))
  expect_equal(getAnnotation(exampleTumourNormal, sampleAnnotation = "tumour", useGroupMeans = FALSE) %>% dim(), c(10,1))
  expect_equal(getAnnotation(exampleTumourNormal, sampleAnnotation = c("tumour","type"),  useGroupMeans = FALSE)  %>% dim(), c(10,2))
  expect_equal(getAnnotation(exampleTumourNormal, sampleAnnotation = c(tumour,type), useGroupMeans = FALSE)  %>% dim(), c(10,2))

  expect_equal(getAnnotation(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = tumour, useGroupMeans = TRUE)  %>% dim(), c(2,1))
  expect_equal(getAnnotation(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = "tumour", useGroupMeans = TRUE)  %>% dim(), c(2,1))

  expect_error(getAnnotation(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = c("tumour","type"), useGroupMeans = TRUE)  %>% dim())
  expect_error(getAnnotation(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = c(tumour,type), useGroupMeans = TRUE)  %>% dim())

  expect_no_error(plotCorrelationMatrix(exampleTumourNormal))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal, sampleAnnotation = tumour, useGroupMeans = FALSE))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal, sampleAnnotation = "tumour", useGroupMeans = FALSE))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal, sampleAnnotation = c("tumour","type"),  useGroupMeans = FALSE))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal, sampleAnnotation = c(tumour,type), useGroupMeans = FALSE))

  expect_no_error(plotCorrelationMatrix(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = tumour, useGroupMeans = TRUE))
  expect_no_error(plotCorrelationMatrix(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = "tumour", useGroupMeans = TRUE))

  expect_error(plotCorrelationMatrix(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = c("tumour","type"), useGroupMeans = TRUE))
  expect_error(plotCorrelationMatrix(exampleTumourNormal %>% mutate(group = tumour), sampleAnnotation = c(tumour,type), useGroupMeans = TRUE))

  regions <- getRegions(exampleTumourNormal)[1:10]
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal, regionsToOverlap = regions, sampleAnnotation = tumour, useGroupMeans = FALSE))
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal, regionsToOverlap = regions, sampleAnnotation = "tumour", useGroupMeans = FALSE))
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal, regionsToOverlap = regions, sampleAnnotation = c("tumour","type"),  useGroupMeans = FALSE))
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal, regionsToOverlap = regions, sampleAnnotation = c(tumour,type), useGroupMeans = FALSE))

  expect_no_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = tumour, useGroupMeans = TRUE))
  expect_no_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = "tumour", useGroupMeans = TRUE))

  expect_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = c("tumour","type"), useGroupMeans = TRUE))
  expect_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = c(tumour,type), useGroupMeans = TRUE))

  expect_no_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = c(tumour,type),
               annotationColors = list(tumour = c("Normal" = "blue", "Tumour" = "firebrick4"))))

  expect_no_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour), regionsToOverlap = regions, sampleAnnotation = c(tumour,type,age),
                                     annotationColors = list(tumour = c("Normal" = "blue", "Tumour" = "firebrick4"))))

  expect_no_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour),
                                     regionsToOverlap = regions %>% mutate(newCol = rnorm(10)),
                                     sampleAnnotation = c(tumour,type,age),
                                     windowAnnotation = c(CpG_density,newCol),
                                     annotationColors = list(tumour = c("Normal" = "blue", "Tumour" = "firebrick4"))))

  expect_no_error(plotRegionsHeatmap(exampleTumourNormal %>% mutate(group = tumour),
                                     regionsToOverlap = regions %>% mutate(newCol = rnorm(10)),
                                     windowAnnotation = c(CpG_density, newCol)))

  expect_error(plotRegionsHeatmap(exampleTumourNormal,
                                  regionsToOverlap = regions %>% 
                                    bind_ranges(regions) %>%
                                    mutate(newCol = rnorm(20)),
                                  windowAnnotation = c(CpG_density, newCol)))
  
}
)

test_that("Testing hg38 related annotation/plotting functions", {

  expect_no_error(plotGeneHeatmap(exampleTumourNormal, gene = "HOXA10", sampleAnnotation = tumour, useGroupMeans = FALSE))
  expect_no_error(plotGeneHeatmap(exampleTumourNormal, gene = "HOXA10", sampleAnnotation = "tumour", useGroupMeans = FALSE))
  expect_no_error(plotGeneHeatmap(exampleTumourNormal, gene = "HOXA10", sampleAnnotation = c("tumour","type"),  useGroupMeans = FALSE))
  expect_no_error(plotGeneHeatmap(exampleTumourNormal, gene = "HOXA10", sampleAnnotation = c(tumour,type), useGroupMeans = FALSE))

  expect_no_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = tumour, useGroupMeans = TRUE))
  expect_no_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = "tumour", useGroupMeans = TRUE))

  expect_no_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = tumour, useGroupMeans = TRUE, showSampleNames = FALSE))

  expect_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = c("tumour","type"), useGroupMeans = TRUE))
  expect_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = c(tumour,type), useGroupMeans = TRUE))

  expect_no_error(plotGeneHeatmap(exampleTumourNormal %>% mutate(group = tumour), gene = "HOXA10", sampleAnnotation = c(tumour,type),
                               annotationColors = list(tumour = c("Normal" = "blue", "Tumour" = "firebrick4"))))

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
                 summariseAcrossWindows(regionsToOverlap = getRegions(.)) %>%
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
  expect_equal(randomSet %>% getBetaTable(useGroupMeans = TRUE) %>% dplyr::select(matches("Tum|Norm")) %>% ncol(), randomSet %>% getSampleGroups2() %>% length())

  expect_equal(randomSet %>% getNRPMTable() %>% nrow(), randomSet %>% getRegions() %>% length())
  expect_equal(randomSet %>% getNRPMTable() %>% dplyr::select(matches("Sim")) %>% ncol(), randomSet %>% getSampleNames() %>% length())
  expect_equal(randomSet %>% getNRPMTable(useGroupMeans = TRUE) %>% dplyr::select(matches("Tum|Norm")) %>% ncol(), randomSet %>% getSampleGroups2() %>% length())

  expect_equal(randomSet %>% getNRPMTable(useGroupMeans = TRUE) %>% dplyr::select(matches("Tum|Norm")) %>% ncol(), randomSet %>% getSampleGroups2() %>% length())

  expect_true("valid_fragments" %in% (randomSet %>% getSampleQCSummary() %>% colnames()))

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
                 annotateWindows(genome = "hg38") %>%
                 summariseDMRsByGene() %>%
                 dim(),c(21,4))

})

test_that("Multiple DMRs", {
  expect_no_error(DMRs <- exampleTumourNormal %>%
    calculateDMRs(variable = "type",
                  contrasts = "all"))
})
