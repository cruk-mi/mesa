test_that("the hg19 chr22 tracks tile chr22 in 100 kb bins (#180)", {

    for (track in list(gc_hg19_chr22_100kb, map_hg19_chr22_100kb)) {
        expect_equal(nrow(track), 513L)
        expect_equal(unique(track$chr), "chr22")
        expect_equal(track$start, seq(1L, by = 100000L, length.out = 513L))
        expect_equal(track$end - track$start + 1L, rep(100000L, 513L))
    }

    # the first 16 Mb of chr22 are N, which HMMcopy marks as -1
    expect_true(all(gc_hg19_chr22_100kb$gc[1:10] == -1))
    expect_true(all(gc_hg19_chr22_100kb$gc[gc_hg19_chr22_100kb$gc >= 0] > 0.2))
    expect_true(all(
        map_hg19_chr22_100kb$map >= 0 & map_hg19_chr22_100kb$map <= 1
    ))
})

test_that("addHMMcopyCNV runs on a MEDIPSData chr22 BAM (#180)", {

    skip_if_not_installed("MEDIPSData")
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

    bamFile <- system.file("extdata", "hESCs.Input.chr22.bam",
        package = "MEDIPSData", mustWork = TRUE)

    qs <- data.frame(
        sample_name = "hESC", file_name = bamFile,
        input_file = bamFile, group = "Input"
    ) %>%
        qsea::createQseaSet("BSgenome.Hsapiens.UCSC.hg19",
            chr.select = "chr22", window_size = 5000) %>%
        addHMMcopyCNV(
            windowSize = 1e5, fragmentLength = 300,
            hmmCopyGC = gc_hg19_chr22_100kb,
            hmmCopyMap = map_hg19_chr22_100kb,
            parallel = FALSE
        )

    cnv <- qsea::getCNV(qs)
    expect_length(cnv, 513L)
    expect_equal(colnames(GenomicRanges::mcols(cnv)), "hESC")
    expect_false(anyNA(GenomicRanges::mcols(cnv)$hESC))
    # one copy-number state for the whole of chr22 in this normal sample
    expect_equal(unique(round(GenomicRanges::mcols(cnv)$hESC, 2)), -0.01)

    # each window got the GC and mappability of its own bin
    inputCounts <- qs@libraries$input_counts
    expect_equal(inputCounts$start, gc_hg19_chr22_100kb$start)
    expect_equal(inputCounts$gc, gc_hg19_chr22_100kb$gc)
    expect_equal(inputCounts$map, map_hg19_chr22_100kb$map)
    expect_equal(rownames(qs@libraries$input_file), "hESC")
})

test_that("addHMMcopyCNV stops when the window size does not match (#180)", {

    skip_if_not_installed("MEDIPSData")
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

    bamFile <- system.file("extdata", "hESCs.Input.chr22.bam",
        package = "MEDIPSData", mustWork = TRUE)

    qs <- data.frame(
        sample_name = "hESC", file_name = bamFile,
        input_file = bamFile, group = "Input"
    ) %>%
        qsea::createQseaSet("BSgenome.Hsapiens.UCSC.hg19",
            chr.select = "chr22", window_size = 5000)

    # each 200 kb window overlaps two of the 100 kb bins
    expect_error(
        addHMMcopyCNV(
            qs, windowSize = 2e5, fragmentLength = 300,
            hmmCopyGC = gc_hg19_chr22_100kb,
            hmmCopyMap = map_hg19_chr22_100kb,
            parallel = FALSE
        ),
        "overlap with multiple windows"
    )
})

test_that("removeCNV() zeroes the copy number but keeps the windows", {

    flat <- removeCNV(exampleTumourNormal)
    cnvValues <- as.matrix(GenomicRanges::mcols(qsea::getCNV(flat)))
    original <- qsea::getCNV(exampleTumourNormal)

    expect_length(qsea::getCNV(flat), length(original))
    expect_equal(colnames(cnvValues), qsea::getSampleNames(exampleTumourNormal))
    expect_true(all(cnvValues == 0))
    # the original is untouched
    expect_true(any(as.matrix(GenomicRanges::mcols(original)) != 0))
})
