test_that("Making a simple qseaSet", {

  skip_on_ci()

  if(!rlang::is_installed("MEDIPSData")){
    skip("MEDIPSData Not installed")
  }

  sampleTable <- data.frame(sample_name = c("Normal1","Tumour1"),
                            group = c("Normal1","Tumour1"),
                            file_name = c(system.file("extdata", "NSCLC_MeDIP_1N_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE),
                                          system.file("extdata", "NSCLC_MeDIP_1T_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE)))

  testSet <- makeQset(sampleTable,
                       BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
                       chrSelect = paste0("chr",20:22),
                       windowSize = 300,
                       CNVwindowSize = 1000000,
                       fragmentType = "Sheared",
                       CNVmethod = "none",
                       coverageMethod = "PairedAndR1s",
                       minMapQual = 10,
                       minInsertSize = 70,
                       maxInsertSize = 1000,
                       minReferenceLength = 30,
                       badRegions = NULL,
                       properPairsOnly = FALSE)

  expect_equal(testSet %>% getSampleNames() %>% length(), 2)
  expect_equal(testSet %>% getRegions() %>% seqinfo() %>% length(), 3)
  expect_equal(testSet %>% getRegions() %>% length(), 541532)

})