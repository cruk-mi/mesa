# The option tests run on synthetic yeast reads (helper-bams.R), which build in
# a few seconds. Two blocks then build from the real MEDIPSData reads, one per
# chromosome naming convention (UCSC and NCBI) and per human BSgenome.

# Yeast CpGs are denser than human ones, so the default maxPatternDensity
# leaves no windows to estimate the background from.
makeYeastQset <- function(
    sampleTable = yeastSampleTable(), CNVmethod = "None", ...
) {
    suppressMessages(
        makeQset(
            sampleTable,
            BSgenome = "BSgenome.Scerevisiae.UCSC.sacCer3",
            chrSelect = "chrI",
            windowSize = 300,
            fragmentType = "Sheared",
            CNVmethod = CNVmethod,
            maxPatternDensity = 100,
            parallel = FALSE,
            ...
        )
    )
}

expectYeastQset <- function(testSet) {
    expect_equal(testSet %>% qsea::getSampleNames() %>% length(), 2)
    expect_equal(testSet %>% qsea::getChrNames(), "chrI")
    expect_equal(testSet %>% qsea::getRegions() %>% width() %>% unique(), 300)
    # sacCer3 chrI is 230218 bp long
    expect_equal(testSet %>% qsea::getRegions() %>% length(), 767)
    expect_equal(
        testSet %>% qsea::getCounts() %>% colSums() %>% unname(),
        testSet %>% qsea::getLibSize()
    )
    # every synthetic pair is a proper pair that passes the filters
    expect_equal(unname(qsea::getLibSize(testSet)), c(5000, 5000))

    expect_true("relH" %in% (
        testSet %>%
            addLibraryInformation() %>%
            qsea::getSampleTable() %>%
            colnames()
    ))
}

test_that("Making a qseaSet from paired reads", {

    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    skip_if_not_installed("Rsamtools")

    testSet <- makeYeastQset(coverageMethod = "PairedAndR1s")

    expectYeastQset(testSet)
    expect_false(testSet@parameters$properPairsOnly)
    expect_equal(testSet@parameters$minMapQual, 10)
})

test_that("Making a qseaSet with qsea coverage method", {

    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    skip_if_not_installed("Rsamtools")

    testSet <- makeYeastQset(coverageMethod = "qseaPaired")

    expectYeastQset(testSet)
})

test_that("Making a qseaSet proper pairs only", {

    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    skip_if_not_installed("Rsamtools")

    testSet <- makeYeastQset(
        coverageMethod = "PairedAndR1s", properPairsOnly = TRUE
    )

    expectYeastQset(testSet)
    expect_true(testSet@parameters$properPairsOnly)
})

test_that("makeQset() rejects GC and map files of the wrong window size", {

    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    skip_if_not_installed("Rsamtools")

    expect_error(
        makeYeastQset(
            CNVmethod = "HMMdefault",
            CNVwindowSize = 1000000,
            hmmCopyGC = mesa::gc_hg38_500kb,
            hmmCopyMap = mesa::map_hg38_500kb
        ),
        "window size should be the same as the CNVwindowSize"
    )
})

test_that("Making a hg19 qseaSet", {

    skip_if_not_installed("MEDIPSData")
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

    bamFiles <- vapply(c("1N", "1T"), function(x) {
        system.file("extdata",
            paste0("NSCLC_MeDIP_", x, "_fst_chr_20_21_22.bam"),
            package = "MEDIPSData", mustWork = TRUE)
    }, character(1))

    sampleTable <- data.frame(sample_name = c("Normal1", "Tumour1"),
        group = c("Normal1", "Tumour1"),
        file_name = unname(bamFiles))

    # chr21 is the smallest chromosome in the BAMs, and holds JAM2
    testSet <- makeQset(sampleTable,
        BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
        chrSelect = "chr21",
        windowSize = 300,
        CNVwindowSize = 1000000,
        fragmentType = "Sheared",
        CNVmethod = "None",
        coverageMethod = "PairedAndR1s",
        minMapQual = 10,
        minInsertSize = 70,
        maxInsertSize = 1000,
        minReferenceLength = 30,
        badRegions = NULL,
        properPairsOnly = FALSE,
        parallel = FALSE)

    chr21Length <- GenomeInfoDb::seqlengths(
        BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
    )[["chr21"]]

    expect_equal(testSet %>% qsea::getSampleNames() %>% length(), 2)
    expect_equal(testSet %>% qsea::getChrNames(), "chr21")
    expect_equal(testSet %>% qsea::getRegions() %>% width() %>% unique(), 300)
    expect_equal(
        testSet %>% qsea::getCounts() %>% colSums() %>% unname(),
        testSet %>% qsea::getLibSize()
    )
    expect_equal(
        testSet %>% qsea::getRegions() %>% length(),
        floor(chr21Length / 300)
    )

    expect_true("relH" %in% (
        testSet %>%
            addLibraryInformation() %>%
            qsea::getSampleTable() %>%
            colnames()
    ))

    localMockBiomart()
    plotGeneHeatmap(testSet, "JAM2")

})

test_that("Making a GRCh38 qseaSet", {

    skip_if_not_installed("MEDIPSData")
    skip_if_not_installed("BSgenome.Hsapiens.NCBI.GRCh38")

    # The MEDIPSData reads are aligned to hg19 with "chr" names. makeQset()
    # drops the prefix to match the NCBI windows, so this checks the GRCh38
    # build path, not the biology: hg19 coordinates on GRCh38 windows.
    # chr20 alone is enough for HMMcopy with 500 kb CNV windows; at 1 Mb it
    # has too few bins for the loess GC correction.
    bamFile <- system.file("extdata", "NSCLC_MeDIP_1N_fst_chr_20_21_22.bam",
        package = "MEDIPSData", mustWork = TRUE)

    sampleTable <- data.frame(sample_name = "Test1",
        group = "Group1",
        file_name = bamFile,
        input_file = bamFile)

    # The BAM also holds chr21 and chr22, which must not reach the CpG
    # enrichment (out-of-bound warnings), and HMMcopy's loess warnings are
    # muffled, so nothing should warn. The GC and map files of the right size
    # (500 kb) must be accepted.
    expect_no_warning(
        GRCh38testSet <- makeQset(sampleTable,
            BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38",
            chrSelect = 20,
            windowSize = 200,
            CNVwindowSize = 500000,
            fragmentType = "cfDNA",
            CNVmethod = "HMMdefault",
            coverageMethod = "PairedAndR1s",
            minMapQual = 10,
            minInsertSize = 70,
            maxInsertSize = 1000,
            minReferenceLength = 30,
            hmmCopyGC = mesa::gc_hg38_500kb,
            hmmCopyMap = mesa::map_hg38_500kb,
            badRegions = NULL,
            properPairsOnly = FALSE)
    )

    chr20Length <- GenomeInfoDb::seqlengths(
        BSgenome::getBSgenome("BSgenome.Hsapiens.NCBI.GRCh38")
    )[["20"]]

    expect_equal(GRCh38testSet %>% getSampleNames() %>% length(), 1)
    expect_equal(GRCh38testSet %>% getChrNames() %>% length(), 1)
    expect_equal(
        GRCh38testSet %>% getRegions() %>% length(),
        floor(chr20Length / 200)
    )

    libTable <- getSampleTable(addLibraryInformation(GRCh38testSet))
    expect_true("relH" %in% colnames(libTable))

    stableTable <- getSampleTable(addHyperStableFraction(GRCh38testSet))
    expect_true("hyperStableEdgar" %in% colnames(stableTable))

    qcSummary <- getSampleQCSummary(GRCh38testSet)
    expect_true("relH" %in% colnames(qcSummary))
    expect_true("valid_fragments" %in% colnames(qcSummary))

})

test_that("calculateCGEnrichment works", {

    skip_if_not_installed("MEDIPSData")
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

    bamFile <- system.file("extdata", "NSCLC_MeDIP_1N_fst_chr_20_21_22.bam",
        package = "MEDIPSData", mustWork = TRUE)

    enr <- calculateCGEnrichment(bamFile,
        BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
        exportPath = NULL,
        extend = 0, shift = 0, uniq = 0,
        chr.select = "chr22", paired = TRUE)

    # nReads is pinned exactly. It is the fragment count this BAM yielded
    # under the MEDIPS implementation, recorded in 9997e94 (2023-01-10) and
    # kept when calculateCGEnrichment() was reimplemented on Rsamtools (#81),
    # so its job is to prove the two implementations agree. A change here is
    # either a real regression or a deliberate change to read filtering
    # (minMapQual, properPairsOnly, chr.select) - not a stale expectation to
    # be updated casually. Re-derive by counting the first mates of properly
    # paired reads on chr22 of the BAM above.
    expect_equal(enr$nReads, 636130)

    # TODO: tolerance is relative, so tolerance = 5 passes for relH anywhere
    # from ~0.5 to ~5 against an expected 3.35. These catch gross breakage
    # only. Tighten when the values can be
    # measured against the Rsamtools implementation.
    expect_equal(enr$relH, 3.353462, tolerance = 5)
    expect_equal(enr$GoGe, 1.631612, tolerance = 5)

})
