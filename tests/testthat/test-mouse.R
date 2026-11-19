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
    # The probe and the connection cost about 5 s, so this is opt-in
    skip_if(Sys.getenv("MESA_LIVE_ENSEMBL") == "",
        "Set MESA_LIVE_ENSEMBL=1 to test against live Ensembl.")
    skip_if_not(ensemblReachable("jul2023.archive.ensembl.org"),
        "Could not reach Ensembl.")

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
    skip_if_not_installed("GO.db")

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
})


test_that("annotateWindows() resolves every way of naming a TxDb alike", {

    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("org.Mm.eg.db")
    skip_if_not_installed("GO.db")

    oldOptions <- options(
        mesa_TxDb = getOption("mesa_TxDb"),
        mesa_annoDb = getOption("mesa_annoDb"),
        mesa_genome = getOption("mesa_genome")
    )
    on.exit(options(oldOptions), add = TRUE)

    mm10 <- getExportedValue(
        "TxDb.Mmusculus.UCSC.mm10.knownGene",
        "TxDb.Mmusculus.UCSC.mm10.knownGene"
    )

    # The annotation itself is run for real above; each extra call costs
    # several seconds, and all that differs here is how the TxDb is named.
    seen <- NULL
    local_mocked_bindings(
        annotatePeak = function(peak, TxDb, annoDb, ...) {
            seen <<- list(TxDb = TxDb, annoDb = annoDb)
            stop("mocked annotatePeak")
        },
        .package = "ChIPseeker"
    )

    resolved <- function(...) {
        seen <<- NULL
        try(annotateWindows(getBetaTable(exampleMouse), ...), silent = TRUE)
        seen
    }

    setMesaTxDb(NULL)
    setMesaAnnoDb(NULL)
    setMesaGenome(NULL)
    byString <- resolved(TxDb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
        annoDb = "org.Mm.eg.db")
    byObject <- resolved(TxDb = mm10, annoDb = "org.Mm.eg.db")

    setMesaTxDb("TxDb.Mmusculus.UCSC.mm10.knownGene")
    setMesaAnnoDb("org.Mm.eg.db")
    byOption <- resolved()

    expect_identical(byString$TxDb, mm10)
    expect_identical(byObject$TxDb, mm10)
    expect_identical(byOption$TxDb, mm10)
    expect_identical(byString$annoDb, "org.Mm.eg.db")
    expect_identical(byOption$annoDb, "org.Mm.eg.db")
})
