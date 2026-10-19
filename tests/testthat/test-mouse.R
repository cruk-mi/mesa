test_that("Mouse plotting", {

    localMockBiomart()

    expect_null(getMart(exampleMouse))
    expect_error(plotGeneHeatmap(exampleMouse, gene = "Fbxl18"),
        "Mmusculus")

    mouseMart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
        host = "https://jul2023.archive.ensembl.org")
    plotGeneHeatmap(exampleMouse, gene = "Fbxl18", mart = mouseMart)
})

test_that("plotGeneHeatmap() prefers an explicit mart to the stored one", {

    localMockBiomart()

    # the mock has no genes for this host, so querying it would fail (#173)
    storedMart <- biomaRt::useMart("ensembl",
        dataset = "mmusculus_gene_ensembl", host = "https://stored.invalid")
    explicitMart <- biomaRt::useMart("ensembl",
        dataset = "mmusculus_gene_ensembl",
        host = "https://jul2023.archive.ensembl.org")

    withStored <- setMart(exampleMouse, storedMart)
    expect_no_error(
        plotGeneHeatmap(withStored, gene = "Fbxl18", mart = explicitMart)
    )
    expect_error(
        plotGeneHeatmap(withStored, gene = "Fbxl18"),
        "Could not retrieve gene information"
    )
})

test_that("Mouse plotting against live Ensembl", {

    skip_on_cran()
    skip_if_offline("jul2023.archive.ensembl.org")

    # exampleMouse stores no Mart, so pass a live one
    liveMart <- tryCatch(
        biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
            host = "https://jul2023.archive.ensembl.org"),
        error = function(e) skip("Could not reach Ensembl.")
    )
    testPlotGeneHeatmap(exampleMouse, gene = "Fbxl18", mart = liveMart)
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
