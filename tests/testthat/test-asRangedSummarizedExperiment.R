test_that("asRangedSummarizedExperiment maps the main slots", {
    qs <- exampleTumourNormal
    se <- asRangedSummarizedExperiment(qs)

    expect_s4_class(se, "RangedSummarizedExperiment")
    expect_equal(ncol(se), length(getSampleNames(qs)))
    expect_equal(nrow(se), length(getRegions(qs)))
    expect_equal(colnames(se), getSampleNames(qs))
    expect_equal(SummarizedExperiment::assayNames(se), "counts")
    expect_equal(
        nrow(SummarizedExperiment::colData(se)),
        nrow(getSampleTable(qs))
    )
    expect_true("parameters" %in% names(S4Vectors::metadata(se)))
})

test_that("library and enrichment information is carried over", {
    qs <- exampleTumourNormal
    se <- asRangedSummarizedExperiment(qs)
    cd <- SummarizedExperiment::colData(se)

    # same columns as addLibraryInformation() puts on the sample table
    libCols <- setdiff(
        colnames(getSampleTable(addLibraryInformation(qs))),
        colnames(getSampleTable(qs))
    )
    expect_true(all(libCols %in% colnames(cd)))
    # bin-level enrichment model stays in metadata
    expect_true(
        all(c("factors", "density", "n") %in%
            names(S4Vectors::metadata(se)$enrichment))
    )
})

test_that("matrix library metrics don't cause an issue (e.g. qsea::addCoverage)", {
    qs <- exampleTumourNormal
    libraries <- getLibraries(qs)
    libraries[["file_name"]] <- as.matrix(libraries[["file_name"]])

    qs <- setLibraries(qs, libraries)

    se <- asRangedSummarizedExperiment(qs)

    expect_true("total_fragments" %in%
                colnames(SummarizedExperiment::colData(se)))
})

test_that("sample names ending in assay suffixes are preserved", {
    for (method in c("counts", "nrpm", "beta")) {
        qs <- exampleTumourNormal %>%
            mutate(newName = paste0(getSampleNames(.), "_", method)) %>%
            renameSamples("newName")

        se <- asRangedSummarizedExperiment(qs, assays = method)

        expect_equal(colnames(se), getSampleNames(qs))
        expect_equal(SummarizedExperiment::assayNames(se), method)
    }
})

test_that("optional assays are added on request", {
    se <- asRangedSummarizedExperiment(
        exampleTumourNormal,
        assays = c("counts", "nrpm", "beta", "cnv")
    )
    expect_setequal(
        SummarizedExperiment::assayNames(se),
        c("counts", "nrpm", "beta", "cnv")
    )
    beta <- SummarizedExperiment::assay(se, "beta")
    expect_true(all(beta >= 0 & beta <= 1, na.rm = TRUE))
})

test_that("CNV is expanded to windows, with NA where uncovered", {
    skip_if(
        length(getCNV(exampleTumourNormal)) == 0,
        "example data has no CNV"
    )
    se <- asRangedSummarizedExperiment(exampleTumourNormal, assays = c("counts", "cnv"))
    cnv <- SummarizedExperiment::assay(se, "cnv")
    expect_equal(dim(cnv), dim(se))
})

test_that("CNV bins stay in metadata when the CNV assay is omitted", {
    cnvBins <- getCNV(exampleTumourNormal)

    se <- asRangedSummarizedExperiment(exampleTumourNormal)

    expect_equal(S4Vectors::metadata(se)$cnvBins, cnvBins)
    expect_false("cnv" %in% SummarizedExperiment::assayNames(se))
})

test_that(".expandCNV picks the largest-overlap bin and leaves gaps NA", {
    windows <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1, 101, 201, 301), width = 100)
    )
    cnv <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1, 151), c(150, 250))
    )
    S4Vectors::mcols(cnv)$s1 <- c(1, 3)

    out <- .expandCNV(cnv, windows, "s1")

    # window 2 (101-200): 50bp in each bin, tie -> first bin; window 3 (201-300)
    # overlaps only the second bin; window 4 (301-400) is uncovered.
    expect_equal(unname(out[, "s1"]), c(1, 1, 3, NA))
})

test_that("non-qseaSet input is rejected", {
    expect_error(asRangedSummarizedExperiment(1:3), "qseaSet")
})