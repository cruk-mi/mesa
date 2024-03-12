test_that("Mouse annotation and plotting", {
  
  expect_true("Mart" %in% class(exampleMouse %>% getMart()))

  expect_no_error(plotGeneHeatmap(exampleMouse, 
                                  gene = "Fbxl18")
                  )
  
  #specify a different Mart:
  expect_no_error(plotGeneHeatmap(exampleMouse, 
                                  gene = "Fbxl18",
                                  mart = biomaRt::useMart('ensembl', dataset='mmusculus_gene_ensembl', host = "https://jul2023.archive.ensembl.org")
                                  ))
  
  exampleMouse2 <- exampleMouse %>%
    setMart(biomaRt::useMart('ensembl', dataset='mmusculus_gene_ensembl', host = "https://jul2023.archive.ensembl.org"))
    
  expect_no_error(plotGeneHeatmap(exampleMouse2, 
                                  gene = "Fbxl18")
                  )
  
  
  #expect an error if global settings not set
  setMesaTxDb(NULL)
  setMesaAnnoDb(NULL)
  setMesaGenome(NULL)
  expect_error(exampleMouse %>% 
                 getBetaTable() %>% 
                 annotateWindows()
               )
  
  # should work if TxDb and annoDb are manually specified
  tab1 <- expect_no_error(exampleMouse %>%
                           getBetaTable() %>%
                           annotateWindows(TxDb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
                                           annoDb = "org.Mm.eg.db"))
  
  expect_true(any(stringr::str_detect(tab1$SYMBOL,"Fbxl18")))
  
  # should work if TxDb is specified using the full object not a string
  tab2 <- expect_no_error(exampleMouse %>%
                           getBetaTable() %>%
                           annotateWindows(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                                           annoDb = "org.Mm.eg.db"))
  
  expect_true(any(stringr::str_detect(tab2$SYMBOL,"Fbxl18")))
  
  # should also work if TxDb and annoDb are specified globally
  setMesaTxDb("TxDb.Mmusculus.UCSC.mm10.knownGene")
  setMesaAnnoDb("org.Mm.eg.db")
  tab3 <- expect_no_error(exampleMouse %>%
                           getBetaTable() %>%
                           annotateWindows())
  
  expect_true(any(stringr::str_detect(tab3$SYMBOL,"Fbxl18")))
  expect_equal(tab1, tab2)
  expect_equal(tab1, tab3)
})
