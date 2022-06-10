context("DMR related tests")

test_that("Calculating DMRs", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000000) %>%
    mutateQset(patient = stringr::str_remove(sample_name,"[TN]$"),
               variableWithOneLevel = "Test",
               experiment = ifelse( stringr::str_detect(sample_name,"[1234]"),"A","B"),
               experimentConfounded = ifelse( stringr::str_detect(sample_name,"[1234]N"),"A","B"),
               )

  expect_error( DMRdata <- randomSet %>%
                 calculateDMRs(variable = "group",
                               contrastsToDo = tibble::tibble(sample1 = "Tumor", sample2 = "Normal")), NA) #expect no error!

  expect_true(DMRdata %>% tibble::has_name(c("Tumor_vs_Normal_log2FC","Tumor_vs_Normal_adjPval","Tumor_vs_Normal_betaDelta")) %>% all())

  expect_error( annotatedData <- DMRdata  %>% annotateData(), NA) #expect no error!

  expect_true(annotatedData %>% tibble::has_name(c("SYMBOL","annotation","geneId","geneChr")) %>% all())

  expect_error(DMRdata <- randomSet %>%
                  calculateDMRs(variable = "group", covariates = "patient",
                                contrastsToDo = tibble::tibble(sample1 = "Tumor", sample2 = "Normal")), NA) #expect no error!
})

test_that("plotting DMRs", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000000) %>%
    mutateQset(patient = stringr::str_remove(sample_name,"[TN]$"),
               variableWithOneLevel = "Test",
               experiment = ifelse( stringr::str_detect(sample_name,"[1234]"),"A","B"),
               experimentConfounded = ifelse( stringr::str_detect(sample_name,"[1234]N"),"A","B"),
    )

  DMRdata <- randomSet %>%
                  calculateDMRs(variable = "group",
                                contrastsToDo = tibble::tibble(sample1 = "Tumor", sample2 = "Normal"))


  expect_error( randomSet %>%
    plotGRangesHeatmap(DMRdata %>% filter(abs(Tumor_vs_Normal_log2FC) > 1) ) , NA) #expect no error

  expect_error( randomSet %>%
                  plotGRangesHeatmap(DMRdata %>% dplyr::filter(abs(Tumor_vs_Normal_log2FC) > 1),
                                     normMethod = "beta") , NA)

  expect_error( randomSet %>%
                  plotGRangesHeatmap(DMRdata %>% dplyr::filter(abs(Tumor_vs_Normal_log2FC) > 1),
                                     clusterRows = TRUE) , NA)
  })

