test_that("makeTransposedTable() has a row per sample, a column per window", {

    tbl <- suppressMessages(
        makeTransposedTable(exampleTumourNormal, "nrpm", tumour, tissue)
    )
    windows <- qsea::getRegions(exampleTumourNormal)

    expect_equal(tbl$sample_name, qsea::getSampleNames(exampleTumourNormal))
    expect_equal(
        colnames(tbl)[1:3], c("sample_name", "tumour", "tissue")
    )
    expect_equal(ncol(tbl), 3 + length(windows))
    expect_equal(
        colnames(tbl)[4],
        paste0("chr", GenomicRanges::seqnames(windows)[1], "_",
            GenomicRanges::start(windows)[1], "_",
            GenomicRanges::end(windows)[1])
    )
    expect_equal(tbl$tumour[tbl$sample_name == "Lung2_N"], "Normal")
})

test_that("countWindowsAboveCutoff() counts windows per sample", {

    windows <- qsea::getRegions(exampleTumourNormal)[1:200]
    counts <- suppressMessages(
        countWindowsAboveCutoff(exampleTumourNormal, windows, cutoff = 1)
    )

    expect_equal(nrow(counts), 10)
    expect_true(all(
        c("sample_name", "numOverCutoff", "totalWindowsUsed", "tumour")
        %in% colnames(counts)
    ))
    expect_true(all(counts$totalWindowsUsed == 200))
    # the cutoff is inclusive, as in the nrpm table of the same windows
    nrpm <- suppressMessages(qsea::makeTable(
        exampleTumourNormal, norm_methods = "nrpm", samples = counts$sample_name
    ))[1:200, ]
    expect_equal(
        counts$numOverCutoff,
        unname(colSums(nrpm[, paste0(counts$sample_name, "_nrpm")] >= 1))
    )
    expect_true(all(counts$numOverCutoff > 0))

    # a huge cutoff leaves nothing, and samples can be chosen
    none <- suppressMessages(countWindowsAboveCutoff(
        exampleTumourNormal, windows, samples = c("Colon1_T", "Lung1_N"),
        cutoff = 1e6
    ))
    expect_equal(none$sample_name, c("Colon1_T", "Lung1_N"))
    expect_equal(none$numOverCutoff, c(0, 0))
})

test_that("subsetWindowsOverBackground() splits windows by background", {

    removeWindows <- function(...) {
        capture.output(out <- suppressMessages(
            subsetWindowsOverBackground(exampleTumourNormal, ...)
        ))
        out
    }
    nWindows <- length(qsea::getRegions(exampleTumourNormal))
    nKept <- function(...) length(qsea::getRegions(removeWindows(...)))
    # numWindows sets the background rate: with the 819 windows of the example
    # every window is below background, so assume a genome of 1e6 windows
    samples <- c("Colon1_T", "Colon2_T", "Lung1_T")

    above <- removeWindows(
        keepAbove = TRUE, samples = samples, numWindows = 1e6
    )
    below <- removeWindows(
        keepAbove = FALSE, samples = samples, numWindows = 1e6
    )

    expect_s4_class(above, "qseaSet")
    expect_equal(length(qsea::getRegions(above)), 121)
    expect_equal(length(qsea::getRegions(below)), nWindows - 121)
    # more samples above background, a stricter FDR or a smaller genome keep
    # fewer windows
    kept <- function(...) {
        nKept(keepAbove = TRUE, samples = samples, ...)
    }
    expect_equal(kept(numWindows = 1e6, numAbove = 2), 75)
    expect_equal(kept(numWindows = 1e6, numAbove = 3), 36)
    expect_equal(kept(numWindows = 1e6, FDRthres = 1e-10), 71)
    expect_equal(kept(numWindows = 1e5), 16)
    expect_equal(
        qsea::getSampleNames(above), qsea::getSampleNames(exampleTumourNormal)
    )

    # samples can be given by name or by a regular expression
    expect_message(
        capture.output(subsetWindowsOverBackground(
            exampleTumourNormal, samples = "_T$", FDRthres = 0.5
        )),
        "Filtering out 5 samples containing _T\\$ in the name."
    )
    expect_error(
        subsetWindowsOverBackground(
            exampleTumourNormal, samples = c("Colon1_T", "nope")
        ),
        "Sample nope not present in the qseaSet!"
    )
})

test_that("removeLibraryFactors() sets every library factor to 1", {

    unity <- removeLibraryFactors(exampleTumourNormal)
    expect_true(all(
        unity@libraries$file_name[, "library_factor"] == 1
    ))

    # a library_factor column in the sample table is reset too
    withFactor <- exampleTumourNormal %>%
        mutate(library_factor = seq_along(sample_name)) %>%
        removeLibraryFactors()
    expect_true(all(qsea::getSampleTable(withFactor)$library_factor == 1))
})

test_that("writeBigWigs() writes one bigWig per sample or per group", {

    skip_if_not_installed("rtracklayer")

    folder <- tempfile("bigwigs")
    on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    small <- subsetQset(
        exampleTumourNormal, samplesToKeep = c("Colon1_T", "Colon1_N")
    )

    expect_message(writeBigWigs(small, folder), "Colon1_T")
    expect_setequal(
        list.files(folder), c("Colon1_T_nrpm.bw", "Colon1_N_nrpm.bw")
    )
    imported <- rtracklayer::import(file.path(folder, "Colon1_T_nrpm.bw"))
    expect_gt(length(imported), 0)

    groupFolder <- tempfile("groupBigwigs")
    on.exit(unlink(groupFolder, recursive = TRUE), add = TRUE)
    # (the groups of exampleTumourNormal hold one sample each)
    suppressMessages(writeBigWigs(small, groupFolder, useGroupMeans = TRUE))
    expect_setequal(
        list.files(groupFolder), c("Colon1_T_nrpm.bw", "Colon1_N_nrpm.bw")
    )
})

test_that("getGenomicFeatureDistribution() counts windows per feature", {

    # annotateWindows() is tested on real annotation in test-DMRs.R and
    # test-mouse.R, and takes about 20 s to load it, so a stand-in labels
    # the windows here.
    local_mocked_bindings(
        annotateWindows = function(dataTable, ...) {
            dataTable %>%
                mutate(
                    landscape = rep_len(c("CpG island", "Open sea"), n()),
                    shortAnno = rep_len(
                        c("Promoter", "Intron", "Distal Intergenic"), n()
                    )
                )
        }
    )
    samples <- c("Colon1_T", "Lung1_N")
    small <- subsetQset(exampleTumourNormal, samplesToKeep = samples)

    distribution <- getGenomicFeatureDistribution(small, cutoff = 1)

    labelled <- getDataTable(
        small, normMethod = "nrpm", minEnrichment = 3, addMethodSuffix = TRUE
    ) %>%
        mutate(
            landscape = rep_len(c("CpG island", "Open sea"), n()),
            shortAnno = rep_len(
                c("Promoter", "Intron", "Distal Intergenic"), n()
            )
        )
    islands <- labelled$landscape == "CpG island"

    expect_true(all(
        c("sample_name", "landscape", "shortAnno", "nWindows", "sum",
            "nOverCutoff", "tumour") %in% colnames(distribution)
    ))
    expect_setequal(distribution$sample_name, samples)
    expect_setequal(
        na.omit(distribution$landscape), c("CpG island", "Open sea")
    )
    expect_setequal(
        na.omit(distribution$shortAnno),
        c("Promoter", "Intron", "Distal Intergenic")
    )

    island <- distribution %>%
        filter(landscape == "CpG island", sample_name == "Colon1_T")
    colonIslands <- labelled$Colon1_T_nrpm[islands]
    expect_equal(island$nWindows, sum(!is.na(colonIslands)))
    expect_equal(island$sum, sum(colonIslands, na.rm = TRUE))
    expect_equal(island$nOverCutoff, sum(colonIslands >= 1, na.rm = TRUE))
    promoter <- distribution %>%
        filter(shortAnno == "Promoter", sample_name == "Lung1_N")
    inPromoter <- labelled$shortAnno == "Promoter"
    expect_equal(
        promoter$nWindows, sum(!is.na(labelled$Lung1_N_nrpm[inPromoter]))
    )
})
