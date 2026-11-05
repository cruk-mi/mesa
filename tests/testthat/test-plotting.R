test_that("ensemblReachable() is TRUE when the BioMart service answers", {

    local_mocked_bindings(martServiceStatus = function(url, timeout) 200L)

    expect_true(ensemblReachable("jul2022.archive.ensembl.org"))
})

test_that("ensemblReachable() is FALSE on an HTTP error or a timeout", {

    local_mocked_bindings(martServiceStatus = function(url, timeout) 500L)
    expect_false(ensemblReachable("jul2022.archive.ensembl.org"))

    # a timeout is an error from curlGetHeaders()
    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        stop("libcurl error code 28:\n\tOperation timed out")
    })
    expect_false(ensemblReachable("jul2022.archive.ensembl.org"))
})

test_that("ensemblReachable() probes the BioMart service, not the site root", {

    seen <- NULL
    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        seen <<- list(url = url, timeout = timeout)
        200L
    })

    # the jul2023 root answered while its BioMart service hung (#183)
    service <- "/biomart/martservice?type=registry"
    ensemblReachable("jul2023.archive.ensembl.org", timeout = 2)
    expect_equal(seen$url,
        paste0("https://jul2023.archive.ensembl.org", service))
    expect_equal(seen$timeout, 2)

    ensemblReachable("https://jul2022.archive.ensembl.org")
    expect_equal(seen$url,
        paste0("https://jul2022.archive.ensembl.org", service))
    expect_equal(seen$timeout, 5)
})

test_that("ensemblReachable() returns FALSE at once when the probe errors", {

    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        stop("Connection timed out")
    })

    elapsed <- system.time(reachable <- ensemblReachable("10.255.255.1"))
    expect_false(reachable)
    expect_lt(elapsed[["elapsed"]], 1)
})

test_that("ensemblReachable() rounds the timeout up to whole seconds", {

    # 0 means a 60 s default to curlGetHeaders()
    sent <- NULL
    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        sent <<- timeout
        200L
    })

    ensemblReachable("jul2022.archive.ensembl.org", timeout = 0.5)
    expect_equal(sent, 1)

    ensemblReachable("jul2022.archive.ensembl.org", timeout = 2.2)
    expect_equal(sent, 3)
})

test_that("ensemblReachable() rejects a timeout that is not positive", {

    local_mocked_bindings(martServiceStatus = function(url, timeout) 200L)

    for (bad in list(0, -1, NA_real_, "5", c(1, 2))) {
        expect_error(
            ensemblReachable("jul2022.archive.ensembl.org", timeout = bad),
            "positive number"
        )
    }
})

test_that("plotGeneHeatmap() fails fast when the default archive hangs", {

    # purrr::possibly() swallows errors from useMart(), so record the call
    called <- FALSE
    local_mocked_bindings(ensemblReachable = function(host, ...) FALSE)
    local_mocked_bindings(
        useMart = function(...) {
            called <<- TRUE
            stop("useMart() must not be called")
        },
        .package = "biomaRt"
    )

    expect_error(
        plotGeneHeatmap(exampleTumourNormal, gene = "HOXA9"),
        "did not answer within 5 s"
    )
    expect_false(called)
})

test_that("plotGeneHeatmap() waits for mesa_ensemblTimeout seconds", {

    oldOptions <- options(mesa_ensemblTimeout = 12)
    on.exit(options(oldOptions), add = TRUE)

    waited <- NULL
    local_mocked_bindings(ensemblReachable = function(host, timeout) {
        waited <<- timeout
        FALSE
    })

    expect_error(
        plotGeneHeatmap(exampleTumourNormal, gene = "HOXA9"),
        "did not answer within 12 s"
    )
    expect_equal(waited, 12)
})

test_that("ensemblReachable() probes a mart's host without doubling the path", {

    seen <- NULL
    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        seen <<- url
        200L
    })

    # Mart@host already ends in the BioMart service path
    ensemblReachable("https://x.org:443/biomart/martservice")
    expect_equal(seen, "https://x.org:443/biomart/martservice?type=registry")
})

test_that("plotGeneHeatmap() probes an explicit mart before getBM()", {

    martHost <- "https://jul2023.archive.ensembl.org:443/biomart/martservice"
    explicitMart <- methods::new("Mart", biomart = "ensembl",
        dataset = "mmusculus_gene_ensembl", host = martHost)

    # purrr::possibly() swallows errors from getBM(), so record the call
    probed <- NULL
    queried <- FALSE
    local_mocked_bindings(ensemblReachable = function(host, timeout) {
        probed <<- host
        FALSE
    })
    local_mocked_bindings(
        getBM = function(...) {
            queried <<- TRUE
            stop("getBM() must not be called")
        },
        .package = "biomaRt"
    )

    expect_error(
        plotGeneHeatmap(exampleMouse, gene = "Fbxl18", mart = explicitMart),
        "did not answer within 5 s"
    )
    expect_equal(probed, martHost)
    expect_false(queried)
})

test_that("plotGeneHeatmap() does not probe a mart with its own http_config", {

    localMockBiomart()

    # getBM() uses http_config (proxy, certificates); the probe cannot
    martHost <- "https://jul2023.archive.ensembl.org:443/biomart/martservice"
    proxyMart <- methods::new("Mart", biomart = "ensembl",
        dataset = "mmusculus_gene_ensembl", host = martHost,
        http_config = list(proxy = "http://proxy.invalid"))

    probed <- FALSE
    local_mocked_bindings(ensemblReachable = function(host, timeout) {
        probed <<- TRUE
        FALSE
    })

    expect_no_error(
        plotGeneHeatmap(exampleMouse, gene = "Fbxl18", mart = proxyMart)
    )
    expect_false(probed)
})

test_that("plotDMRUpset() plots overlaps between contrasts", {

    skip_if_not_installed("UpSetR")

    dmrs <- tibble::tibble(
        A_vs_B_adjPval = c(0.001, 0.01, 0.5, 0.01),
        A_vs_C_adjPval = c(0.001, 0.5, 0.01, 0.01),
        B_vs_C_adjPval = c(0.5, 0.5, 0.5, 0.01)
    )

    # UpSetR still uses ggplot2 features deprecated since 3.4
    plot <- suppressWarnings(plotDMRUpset(dmrs))
    expect_s3_class(plot, "upset")
    expect_equal(sort(colnames(plot$New_data)),
        c("A_vs_B", "A_vs_C", "B_vs_C"))
    # four windows, three significant in A_vs_B, three in A_vs_C
    expect_equal(sum(plot$New_data$A_vs_B), 3)
    expect_equal(sum(plot$New_data$A_vs_C), 3)

    # removeVS drops "_vs_..." from the names, which must stay unique
    distinct <- dmrs
    colnames(distinct) <- c(
        "A_vs_B_adjPval", "C_vs_B_adjPval", "D_vs_B_adjPval"
    )
    named <- suppressWarnings(plotDMRUpset(
        distinct, string = "^[AC]$", removeVS = TRUE, minAdjPval = 0.005
    ))
    expect_equal(sort(colnames(named$New_data)), c("A", "C"))
    # only windows below minAdjPval count
    expect_equal(sum(named$New_data$A), 1)
})

test_that("plotDMRUpset() stops without UpSetR or with one contrast", {

    skip_if_not_installed("UpSetR")

    dmrs <- tibble::tibble(
        A_vs_B_adjPval = c(0.001, 0.01), A_vs_C_adjPval = c(0.5, 0.01)
    )
    expect_error(
        plotDMRUpset(dmrs, string = "A_vs_B"),
        "Only 1 column remaining, nothing to plot!"
    )

    local_mocked_bindings(requireNamespace = function(...) FALSE,
        .package = "base")
    expect_error(plotDMRUpset(dmrs), "must be installed")
})

test_that("plotGenomicFeatureDistribution() plots feature fractions", {

    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    oldGenome <- getOption("mesa_genome")
    on.exit(options(mesa_genome = oldGenome), add = TRUE)

    # Annotating peaks with ChIPseeker takes about 10 s to load the
    # annotation, so a stand-in labels the windows with three features.
    # getGenomicFeatureDistribution() tests the real annotation.
    methods::setClass("mockPeakAnno", methods::representation(anno = "GRanges"))
    on.exit(methods::removeClass("mockPeakAnno"), add = TRUE)
    features <- c("Promoter (<=1kb)", "Intron (ENSG1, intron 1 of 3)",
        "Distal Intergenic")
    seen <- list()
    local_mocked_bindings(
        annotatePeak = function(peak, TxDb, annoDb, ...) {
            seen <<- list(TxDb = TxDb, annoDb = annoDb)
            peak$annotation <- rep_len(features, length(peak))
            methods::new("mockPeakAnno", anno = peak)
        },
        .package = "ChIPseeker"
    )
    small <- subsetQset(
        exampleTumourNormal, samplesToKeep = c("Colon1_T", "Colon1_N")
    )
    nrpm <- getNRPMTable(small)
    expectedBars <- vapply(
        c("Colon1_T", "Colon1_N"),
        function(x) sum(nrpm[[x]] > 1), numeric(1)
    )

    setMesaGenome(NULL)
    plot <- plotGenomicFeatureDistribution(
        small, cutoff = 1, barType = "fill", genome = "hg38"
    )

    expect_s3_class(plot, "ggplot")
    expect_setequal(unique(as.character(plot$data$sample)),
        c("Colon1_T", "Colon1_N"))
    # "Intron (ENSG1, intron 1 of 3)" is shortened to its feature
    expect_setequal(unique(plot$data$feature),
        c("Promoter (<=1kb)", "Intron", "Distal Intergenic"))
    expect_equal(
        vapply(c("Colon1_T", "Colon1_N"),
            function(x) sum(plot$data$value[plot$data$sample == x]),
            numeric(1)),
        expectedBars
    )
    expect_equal(plot$layers[[1]]$position$fill, TRUE)
    expect_s4_class(seen$TxDb, "TxDb")
    expect_equal(seen$annoDb, "org.Hs.eg.db")

    # the global genome, then the hg38 fallback, give the same annotation
    setMesaGenome("hg38")
    viaOption <- plotGenomicFeatureDistribution(small, cutoff = 1)
    setMesaGenome(NULL)
    fallback <- plotGenomicFeatureDistribution(small, cutoff = 1)
    expect_equal(viaOption$data, plot$data)
    expect_equal(fallback$data, plot$data)
    expect_equal(seen$annoDb, "org.Hs.eg.db")

    expect_error(
        plotGenomicFeatureDistribution(small, genome = "hg5"),
        "Unsupported genome: hg5"
    )
})
