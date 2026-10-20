test_that("Mouse plotting", {

    localMockBiomart()

    expect_true("Mart" %in% class(exampleMouse %>% getMart()))

    plotGeneHeatmap(exampleMouse, gene = "Fbxl18")
})

test_that("Mouse plotting against live Ensembl", {

    skip_on_cran()
    skip_if_offline("jul2023.archive.ensembl.org")

    # exampleMouse stores a Mart for a retired archive, so set a live one
    liveMart <- tryCatch(
        biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
            host = "https://jul2023.archive.ensembl.org"),
        error = function(e) skip("Could not reach Ensembl.")
    )
    testPlotGeneHeatmap(exampleMouse %>% setMart(liveMart), gene = "Fbxl18")
})

test_that("Mouse annotation", {

    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("org.Mm.eg.db")

    oldOptions <- options(
        mesa_TxDb = getOption("mesa_TxDb"),
        mesa_annoDb = getOption("mesa_annoDb"),
        mesa_genome = getOption("mesa_genome")
    )
    on.exit(options(oldOptions), add = TRUE)

    # expect an error if global settings not set
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

    expect_true(any(stringr::str_detect(tab1$SYMBOL, "Fbxl18")))

    # should work if TxDb is specified using the full object not a string
    tab2 <- expect_no_error(exampleMouse %>%
        getBetaTable() %>%
        annotateWindows(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
            annoDb = "org.Mm.eg.db"))

    expect_true(any(stringr::str_detect(tab2$SYMBOL, "Fbxl18")))

    # should also work if TxDb and annoDb are specified globally
    setMesaTxDb("TxDb.Mmusculus.UCSC.mm10.knownGene")
    setMesaAnnoDb("org.Mm.eg.db")
    tab3 <- expect_no_error(exampleMouse %>%
        getBetaTable() %>%
        annotateWindows())

    expect_true(any(stringr::str_detect(tab3$SYMBOL, "Fbxl18")))
    expect_equal(tab1, tab2)
    expect_equal(tab1, tab3)
})
