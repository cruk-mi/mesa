test_that("Calculating DMRs", {

  BiocParallel::SerialParam()

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000000) %>%
    mutate(patient = stringr::str_remove(sample_name,"[TN]$"),
               variableWithOneLevel = "Test",
               experiment = ifelse( stringr::str_detect(sample_name,"[1234]"),"A","B"),
               experimentConfounded = ifelse( stringr::str_detect(sample_name,"[1234]N"),"A","B"),
               )

  expect_no_error( DMRdata <- randomSet %>%
                 calculateDMRs(variable = "group",
                               contrasts = tibble::tibble(sample1 = "Tumor", sample2 = "Normal")))

  expect_true(DMRdata %>% tibble::has_name(c("Tumor_vs_Normal_log2FC","Tumor_vs_Normal_adjPval","Tumor_vs_Normal_deltaBeta")) %>% all())

  expect_no_error(annotatedData <- DMRdata  %>% annotateWindows()) #expect no error!

  expect_true(annotatedData %>% tibble::has_name(c("SYMBOL","annotation","geneId","geneChr")) %>% all())

  expect_no_error(DMRdata <- randomSet %>%
                  calculateDMRs(variable = "group", covariates = "patient",
                                contrasts = tibble::tibble(group1 = "Tumor", group2 = "Normal")))

  expect_no_error(summariseDefault <- randomSet %>% summariseAcrossWindows(DMRdata))

  expect_true("group" %in% colnames(summariseDefault))
  expect_true("beta_mean" %in% colnames(summariseDefault))
  expect_true("nrpm_mean" %in% colnames(summariseDefault))

  expect_no_error(summarise2 <- randomSet %>% summariseAcrossWindows(DMRdata, fn = median, normMethod = "counts", addSampleTable = FALSE))

  expect_false("group" %in% colnames(summarise2))
  expect_true("counts_median" %in% colnames(summarise2))
  expect_false("nrpm_mean" %in% colnames(summarise2))

  expect_no_error(setWithSummary <- randomSet %>%
                    addSummaryAcrossWindows(DMRdata, fn = median, normMethod = "beta") %>%
                    addSummaryAcrossWindows(DMRdata %>% filter(Tumor_vs_Normal_log2FC >= 2), fn = max, normMethod = "nrpm", suffix = "highest_only"))

  expect_true(setWithSummary %>% is.qseaSet())
  expect_true("beta_median" %in% colnames(getSampleTable(setWithSummary)))
  expect_true("nrpm_max_highest_only" %in% colnames(getSampleTable(setWithSummary)))

  expect_equal( exampleTumourNormal %>%
                  mutate(group = stringr::str_remove(group,"\\d")) %>%
                  calculateDMRs(variable = "tumour",
                                contrasts = "Tumour_vs_Normal",
                                keepContrastMeans = TRUE,
                                keepData = FALSE,
                                keepGroupMeans = FALSE) %>%
                  ncol(), 11)

  expect_equal( exampleTumourNormal %>%
                  mutate(group = stringr::str_remove(group,"\\d")) %>%
                     calculateDMRs(variable = "tumour",
                                   contrasts = "Tumour_vs_Normal",
                                   keepContrastMeans = TRUE,
                                   keepData = TRUE,
                                   keepGroupMeans = FALSE) %>%
                     ncol(), 31)

  expect_equal( exampleTumourNormal %>%
                  mutate(group = stringr::str_remove(group,"\\d")) %>%
                     calculateDMRs(variable = "tumour",
                                   contrasts = "Tumour_vs_Normal",
                                   keepContrastMeans = TRUE,
                                   keepData = FALSE,
                                   keepGroupMeans = TRUE) %>%
                     ncol(), 19)

  expect_equal(exampleTumourNormal %>%
                 mutate(group = stringr::str_remove(group,"\\d")) %>%
                     calculateDMRs(variable = "tumour",
                                   contrasts = "Tumour_vs_Normal",
                                   keepContrastMeans = TRUE,
                                   keepData = TRUE,
                                   keepGroupMeans = TRUE) %>%
                    ncol(), 39)

  expect_equal( exampleTumourNormal %>%
                  mutate(group = stringr::str_remove(group,"\\d")) %>%
                  calculateDMRs(variable = "tumour",
                                contrasts = "Tumour_vs_Normal",
                                keepContrastMeans = FALSE,
                                keepData = FALSE,
                                keepGroupMeans = FALSE) %>%
                  ncol(), 7)

  expect_equal( exampleTumourNormal %>%
                  mutate(group = stringr::str_remove(group,"\\d")) %>%
                  calculateDMRs(variable = "tumour",
                                contrasts = "Tumour_vs_Normal",
                                keepContrastMeans = FALSE,
                                keepData = TRUE,
                                keepGroupMeans = FALSE) %>%
                  ncol(), 27)

  expect_equal( exampleTumourNormal %>%
                  mutate(group = stringr::str_remove(group,"\\d")) %>%
                  calculateDMRs(variable = "tumour",
                                contrasts = "Tumour_vs_Normal",
                                keepContrastMeans = FALSE,
                                keepData = FALSE,
                                keepGroupMeans = TRUE) %>%
                  ncol(), 15)

  expect_equal(exampleTumourNormal %>%
                 mutate(group = stringr::str_remove(group,"\\d")) %>%
                 calculateDMRs(variable = "tumour",
                               contrasts = "Tumour_vs_Normal",
                               keepContrastMeans = FALSE,
                               keepData = TRUE,
                               keepGroupMeans = TRUE) %>%
                 ncol(), 35)

  expect_true(DMRdata %>% tibble::has_name(c("Tumor_vs_Normal_log2FC","Tumor_vs_Normal_adjPval","Tumor_vs_Normal_deltaBeta")) %>% all())

  testWithFilter <- exampleTumourNormal %>%
    filter(type %in% c("LUAD", "NormalLung")) %>%
    calculateDMRs(variable = "type", 
                  contrasts = "LUAD_vs_NormalLung")
  
  testNoFilter <- exampleTumourNormal %>%
    calculateDMRs(variable = "type", 
                  contrasts = "LUAD_vs_NormalLung")
  
  expect_equal(testWithFilter %>% dplyr::select(seqnames, start, end, matches("_vs_")), 
               testNoFilter %>% dplyr::select(seqnames, start, end, matches("_vs_")))
  
  testNoFilterShared <- exampleTumourNormal %>%
    filter(type %in% c("LUAD", "NormalLung")) %>%
    calculateDMRs(variable = "type", 
                  contrasts = "LUAD_vs_NormalLung",
                  calcDispersionAll = TRUE
                  )  
  
  expect_not_equal(testNoFilter %>% nrow(), 
               testNoFilterShared %>% nrow())
  
})

test_that("plotting DMRs", {

  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000) %>%
    mutate(patient = stringr::str_remove(sample_name,"[TN]$"),
           variableWithOneLevel = "Test",
           experiment = ifelse( stringr::str_detect(sample_name,"[1234]"),"A","B"),
           experimentConfounded = ifelse( stringr::str_detect(sample_name,"[1234]N"),"A","B"),
           numeric1 = rnorm(16),
           numeric2 = runif(16),
    )

  DMRdata <- randomSet %>%
                  calculateDMRs(variable = "group",
                                contrasts = tibble::tibble(sample1 = "Tumor", sample2 = "Normal"))


  expect_no_error( randomSet %>%
    plotRegionsHeatmap(DMRdata %>% filter(abs(Tumor_vs_Normal_log2FC) > 1) )) 

  expect_no_error( randomSet %>%
                  plotRegionsHeatmap(DMRdata %>% dplyr::filter(abs(Tumor_vs_Normal_log2FC) > 1),
                                     normMethod = "beta"))

  expect_no_error( randomSet %>%
                  plotRegionsHeatmap(DMRdata %>% dplyr::filter(abs(Tumor_vs_Normal_log2FC) > 1),
                                     clusterRows = TRUE))
  
  expect_no_error(randomSet %>%
                     plotRegionsHeatmap(DMRdata %>% filter(abs(Tumor_vs_Normal_log2FC) > 1),
                                        sampleAnnotation = "experiment",
                                        windowAnnotation = "Tumor_vs_Normal_log2FC")) 
  
  expect_no_error(randomSet %>%
                     plotRegionsHeatmap(DMRdata %>% filter(abs(Tumor_vs_Normal_log2FC) > 1) %>% arrange(CpG_density),
                                        sampleAnnotation = "experiment",
                                        windowAnnotation = "Tumor_vs_Normal_log2FC")) 
  
  expect_no_error(randomSet %>%
                    plotRegionsHeatmap(DMRdata %>% filter(abs(Tumor_vs_Normal_log2FC) > 1),
                                       clusterRows = FALSE,
                                       sampleAnnotation = "experiment",
                                       windowAnnotation = c("CpG_density","Tumor_vs_Normal_log2FC")))
  
  expect_no_error( randomSet %>% 
                  plotRegionsHeatmap(DMRdata %>% dplyr::filter(abs(Tumor_vs_Normal_log2FC) > 1),
                                     sampleAnnotation = c("group", "experiment"),
                                     annotationColors = list(group = c("Tumor" = "blue", "Normal" = "red"))))

  expect_no_error( randomSet %>% 
                  plotRegionsHeatmap(DMRdata %>% dplyr::filter(abs(Tumor_vs_Normal_log2FC) > 1),
                                     sampleAnnotation = c("group", "experiment"),
                                     annotationColors = list(group = c("Tumor" = "blue", "Normal" = "red"), experiment = c("A" = "green", "B" = "orange"))))

  })
