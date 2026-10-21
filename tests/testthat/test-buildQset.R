# Fast tests for the read-counting and CNV helpers behind makeQset(). They use
# a tiny synthetic BAM on the small yeast genome, so unlike the long checks in
# test-makeQset.R they need neither MEDIPSData nor a human BSgenome.

# A coordinate-sorted BAM with proper pairs on chrI and on chrII. The chrII
# header is longer than sacCer3's chrII (813184), and the chrII pair lies past
# its end, as reads from another genome build can.
writeYeastBam <- function(dir) {

    seq <- function(n) strrep("A", n)
    qual <- function(n) strrep("I", n)

    rec <- function(name, flag, chr, pos, rnext, pnext, tlen) {
        paste(
            name, flag, chr, pos, 60, "50M", rnext, pnext, tlen,
            seq(50), qual(50),
            sep = "\t"
        )
    }

    lines <- c(
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chrI\tLN:230218",
        "@SQ\tSN:chrII\tLN:900000",
        # two pairs, midpoints at 10100 (window 1) and 10600 (window 3)
        rec("pe1", 99, "chrI", 10000, "=", 10150, 200),
        rec("pe1", 147, "chrI", 10150, "=", 10000, -200),
        rec("pe2", 99, "chrI", 10500, "=", 10650, 200),
        rec("pe2", 147, "chrI", 10650, "=", 10500, -200),
        # a pair beyond the end of sacCer3's chrII
        rec("pe3", 99, "chrII", 850000, "=", 850150, 200),
        rec("pe3", 147, "chrII", 850150, "=", 850000, -200)
    )

    sam <- file.path(dir, "yeast.sam")
    writeLines(lines, sam)
    Rsamtools::asBam(sam, file.path(dir, "yeast"), overwrite = TRUE)
}

test_that("reads outside the regions' chromosomes are dropped (#122)", {

    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")
    skip_if_not_installed("Rsamtools")

    dir <- tempfile("mesaYeast")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bam <- writeYeastBam(dir)

    # The regions declare both chromosomes (as the windows of a qseaSet do)
    # but only have windows on chrI.
    regions <- GenomicRanges::GRanges(
        "chrI",
        IRanges::IRanges(start = c(10001, 10201, 10401), width = 200),
        seqinfo = GenomeInfoDb::Seqinfo(
            c("chrI", "chrII"), seqlengths = c(230218, 813184)
        )
    )

    expect_no_warning(
        res <- getBamCoveragePairedAndUnpairedR1(
            fileName = bam,
            BSgenome = "BSgenome.Scerevisiae.UCSC.sacCer3",
            regions = regions,
            fragmentLength = 150,
            minMapQual = 10,
            minInsertSize = 50,
            maxInsertSize = 1000,
            minReferenceLength = 30,
            properPairsOnly = FALSE
        )
    )

    expect_equal(res$regionsGR$counts, c(1, 0, 1))
})

test_that("muffleLoessWarnings() drops only HMMcopy's loess warnings", {

    expect_no_warning(
        out <- muffleLoessWarnings({
            warning("span too small.   fewer data values than df")
            warning("pseudoinverse used at 0.36828")
            warning("neighborhood radius 0.0057365")
            warning("reciprocal condition number  0")
            warning("There are other near singularities as well. 0.0036292")
            warning("k-d tree limited by memory. ncmax= 200")
            warning("zero-width neighborhood. make span bigger")
            "value"
        })
    )
    expect_equal(out, "value")

    expect_warning(
        muffleLoessWarnings({
            warning("span too small")
            warning("something else")
        }),
        "something else"
    )
})

test_that("runHMMCopy() segments a table without loess warnings", {

    skip_if_not_installed("HMMcopy")

    set.seed(1)
    n <- 400L
    w <- 5e4
    cnv <- tibble::tibble(
        chr = factor(
            rep(c("chr1", "chr2"), each = n / 2),
            levels = c("chr1", "chr2")
        ),
        start = rep(seq(1, by = w, length.out = n / 2), 2),
        end = start + w - 1L,
        width = w,
        strand = "*",
        gc = stats::runif(n, 0.05, 0.95),
        map = stats::runif(n, 0.60, 0.98),
        sampleA = stats::rpois(n, rep(c(500, 1000, 500), c(300, 60, 40)))
    )

    expect_no_warning(out <- runHMMCopy(cnv, colname = "sampleA"))
    expect_s4_class(out, "GRanges")
    expect_true("sampleA" %in% colnames(GenomicRanges::mcols(out)))
})
