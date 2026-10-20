test_that("Making a hg19 qseaSet", {

    # skip check unless options(run_long_checks = TRUE)
    skip_long_checks()

    if (!rlang::is_installed("MEDIPSData")) {
        skip("MEDIPSData Not installed")
    }

    sampleTable <- data.frame(sample_name = c("Normal1", "Tumour1"),
        group = c("Normal1", "Tumour1"),
        file_name = c(system.file("extdata", "NSCLC_MeDIP_1N_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE),
            system.file("extdata", "NSCLC_MeDIP_1T_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE)))

    testSet <- makeQset(sampleTable,
        BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
        chrSelect = paste0("chr", 20:22),
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

    expect_equal(testSet %>% qsea::getSampleNames() %>% length(), 2)
    expect_equal(testSet %>% qsea::getChrNames() %>% length(), 3)
    expect_equal(testSet %>% qsea::getRegions() %>% width() %>% unique(), 300)
    expect_equal(testSet %>% qsea::getCounts() %>% colSums() %>% unname(), testSet %>% qsea::getLibSize())
    expect_equal(testSet %>% qsea::getRegions() %>% length(), 541532)

    expect_true("relH" %in% (testSet %>% addLibraryInformation() %>% qsea::getSampleTable() %>% colnames()))

    localMockBiomart()
    plotGeneHeatmap(testSet, "JAM2")

})

test_that("Making a hg19 qseaSet with qsea coverage method", {

    # skip check unless options(run_long_checks = TRUE)
    skip_long_checks()

    if (!rlang::is_installed("MEDIPSData")) {
        skip("MEDIPSData Not installed")
    }

    sampleTable <- data.frame(sample_name = c("Test1", "Test2"),
        group = c("Test1", "Test2"),
        file_name = c(system.file("extdata", "NSCLC_MeDIP_2N_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE),
            system.file("extdata", "NSCLC_MeDIP_2T_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE)),
        input_file = c(system.file("extdata", "NSCLC_MeDIP_3N_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE),
            system.file("extdata", "NSCLC_MeDIP_3T_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE))
    )

    testSet <- makeQset(sampleTable,
        BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
        chrSelect = paste0("chr", 22),
        windowSize = 300,
        CNVwindowSize = 1000000,
        fragmentType = "Sheared",
        CNVmethod = "None",
        coverageMethod = "qseaPaired",
        minMapQual = 10,
        minInsertSize = 70,
        maxInsertSize = 1000,
        minReferenceLength = 30,
        badRegions = NULL,
        properPairsOnly = FALSE)

    expect_equal(testSet %>% qsea::getSampleNames() %>% length(), 2)
    expect_equal(testSet %>% qsea::getChrNames() %>% length(), 1)
    expect_equal(testSet %>% qsea::getRegions() %>% width() %>% unique(), 300)
    expect_equal(testSet %>% qsea::getCounts() %>% colSums() %>% unname(), testSet %>% qsea::getLibSize())
    expect_equal(testSet %>% qsea::getRegions() %>% length(), 171015)

    expect_true("relH" %in% (testSet %>% addLibraryInformation() %>% qsea::getSampleTable() %>% colnames()))

    localMockBiomart()
    plotGeneHeatmap(testSet, "EWSR1")

})

test_that("Making a GRCh38 qseaSet", {

    # skip check unless options(run_long_checks = TRUE)
    skip_long_checks()
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
        badRegions = NULL,
        properPairsOnly = FALSE)

    expect_equal(GRCh38testSet %>% getSampleNames() %>% length(), 1)
    expect_equal(GRCh38testSet %>% getChrNames() %>% length(), 1)
    expect_equal(GRCh38testSet %>% getRegions() %>% length(), 322220)

    libTable <- getSampleTable(addLibraryInformation(GRCh38testSet))
    expect_true("relH" %in% colnames(libTable))

    stableTable <- getSampleTable(addHyperStableFraction(GRCh38testSet))
    expect_true("hyperStableEdgar" %in% colnames(stableTable))

    qcSummary <- getSampleQCSummary(GRCh38testSet)
    expect_true("relH" %in% colnames(qcSummary))
    expect_true("valid_fragments" %in% colnames(qcSummary))

    # check error is given if wrong size GC/map files given
    expect_error(
        makeQset(sampleTable,
            BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38",
            chrSelect = 22,
            windowSize = 200,
            CNVwindowSize = 1000000,
            fragmentType = "cfDNA",
            CNVmethod = "HMMdefault",
            coverageMethod = "PairedAndR1s",
            hmmCopyGC = mesa::gc_hg38_500kb,
            hmmCopyMap = mesa::map_hg38_500kb,
            badRegions = NULL),
        "window size should be the same as the CNVwindowSize"
    )

    # expect no error if right size GC/map files is given
    expect_no_error(
        makeQset(sampleTable,
            BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38",
            chrSelect = 22,
            windowSize = 200,
            CNVwindowSize = 500000,
            fragmentType = "cfDNA",
            CNVmethod = "None",
            coverageMethod = "PairedAndR1s",
            hmmCopyGC = mesa::gc_hg38_500kb,
            hmmCopyMap = mesa::map_hg38_500kb,
            badRegions = NULL)
    )

})

test_that("Making a qseaSet proper pairs only", {

    # skip check unless options(run_long_checks = TRUE)
    skip_long_checks()
    skip_if_not_installed("MEDIPSData")
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

    bamFile <- system.file("extdata", "NSCLC_MeDIP_1N_fst_chr_20_21_22.bam",
        package = "MEDIPSData", mustWork = TRUE)

    sampleTable <- data.frame(sample_name = "Test1",
        group = "Group1",
        file_name = bamFile,
        input_file = bamFile)

    testSet <- makeQset(sampleTable,
        BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
        chrSelect = "chr22",
        windowSize = 200,
        CNVwindowSize = 1000000,
        fragmentType = "cfDNA",
        CNVmethod = "None",
        coverageMethod = "PairedAndR1s",
        minMapQual = 10,
        minInsertSize = 70,
        maxInsertSize = 1000,
        minReferenceLength = 30,
        badRegions = NULL,
        properPairsOnly = TRUE)

    expect_equal(testSet %>% getSampleNames() %>% length(), 1)
    expect_equal(testSet %>% getChrNames(), "chr22")
    expect_equal(testSet %>% getRegions() %>% length(), 256522)

    libTable <- getSampleTable(addLibraryInformation(testSet))
    expect_true("relH" %in% colnames(libTable))

})

test_that("calculateCGEnrichment works", {
    skip_long_checks()

    if (!rlang::is_installed("MEDIPSData")) {
        skip("MEDIPSData Not installed")
    }

    enr <- calculateCGEnrichment(system.file("extdata", "NSCLC_MeDIP_1N_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE),
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
    # only. Tighten when this block is re-enabled and the values can be
    # measured against the Rsamtools implementation.
    expect_equal(enr$relH, 3.353462, tolerance = 5)
    expect_equal(enr$GoGe, 1.631612, tolerance = 5)

})

test_that("calculateCGEnrichment works (single-end)", {

    # skip check unless options(run_long_checks = TRUE)
    skip_long_checks()

    if (!rlang::is_installed("MEDIPSData")) {
        skip("MEDIPSData Not installed")
    }

    # The fresh-session failure (#81) only manifested when GenomicRanges was
    # not attached to the search path, so the assertion below is only
    # meaningful while that holds. It is NOT order-independent: with
    # skip_long_checks disabled the makeQset blocks above run first and
    # attach GenomicRanges, so skip rather than report a false failure.
    # The same guarantee runs unconditionally in test-readFragments.R,
    # which needs neither MEDIPSData nor a BSgenome.
    skip_if(
        "package:GenomicRanges" %in% search(),
        "GenomicRanges already attached by an earlier test"
    )
    expect_false("package:GenomicRanges" %in% search())

    enr <- calculateCGEnrichment(system.file("extdata", "NSCLC_MeDIP_1N_fst_chr_20_21_22.bam", package = "MEDIPSData", mustWork = TRUE),
        BSgenome = "BSgenome.Hsapiens.UCSC.hg19",
        exportPath = NULL,
        extend = 0, shift = 0, uniq = 0,
        chr.select = "chr22", paired = FALSE)

    expect_true(all(c("file", "relH", "GoGe", "nReads",
        "nReadsWithoutPattern", "n100bpReads",
        "n100bpReadsWithoutPattern") %in% colnames(enr)))
    expect_gt(enr$nReads, 0)
    expect_true(is.numeric(enr$relH) && is.finite(enr$relH))
    expect_true(is.numeric(enr$GoGe) && is.finite(enr$GoGe))

})
