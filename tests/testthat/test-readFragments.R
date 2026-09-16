# Unit tests for the internal BAM import helpers that replaced
# MEDIPS::getGRange() and MEDIPS::getPairedGRange() (#81).
#
# These build a tiny synthetic BAM at run time, so unlike the enrichment tests
# in test-makeQset.R they need neither MEDIPSData nor a BSgenome and run under
# a standard R CMD check.

# A minimal coordinate-sorted SAM exercising every filtering rule the helpers
# inherit from MEDIPS.
#
# chr1 carries the single-end material, chr2 the paired-end material, so each
# helper can be tested with chr.select pointing at its own chromosome.
writeTestBam <- function(dir) {

    seq <- function(n) strrep("A", n)
    qual <- function(n) strrep("I", n)

    rec <- function(name, flag, chr, pos, cigar, len,
        rnext = "*", pnext = 0, tlen = 0) {
        paste(
            name, flag, chr, pos, 60, cigar, rnext, pnext, tlen,
            seq(len), qual(len),
            sep = "\t"
        )
    }

    lines <- c(
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:10000",
        "@SQ\tSN:chr2\tLN:10000",

        # --- chr1: single-end ---
        rec("se_long", 0, "chr1", 1000, "100M", 100),
        rec("se_short", 0, "chr1", 2000, "30M", 30),
        rec("se_dupA", 0, "chr1", 3000, "50M", 50),
        rec("se_dupB", 0, "chr1", 3000, "50M", 50),
        # same coordinates as the pair above but on the minus strand: MEDIPS
        # deduplicated before setting strand to "*", so uniq must keep this.
        rec("se_rev", 16, "chr1", 3000, "50M", 50),
        # secondary alignment: excluded by isSecondaryAlignment = FALSE
        rec("se_secondary", 256, "chr1", 4000, "50M", 50),
        # soft-clipped: KEPT, because mesa passes simpleCigar = FALSE
        rec("se_softclip", 0, "chr1", 5000, "10S40M", 50),

        # --- chr2: paired-end, properly paired ---
        rec("pe1", 99, "chr2", 1000, "50M", 50, "=", 1150, 200),
        rec("pe1", 147, "chr2", 1150, "50M", 50, "=", 1000, -200),
        # secondary pair: excluded
        rec("pe_sec", 355, "chr2", 2000, "50M", 50, "=", 2150, 200),
        rec("pe_sec", 403, "chr2", 2150, "50M", 50, "=", 2000, -200),
        # two fragments with identical spans, for uniq
        rec("pe_dupA", 99, "chr2", 3000, "50M", 50, "=", 3150, 200),
        rec("pe_dupB", 99, "chr2", 3000, "50M", 50, "=", 3150, 200),
        rec("pe_dupA", 147, "chr2", 3150, "50M", 50, "=", 3000, -200),
        rec("pe_dupB", 147, "chr2", 3150, "50M", 50, "=", 3000, -200)
    )

    sam <- file.path(dir, "test.sam")
    writeLines(lines, sam)

    # asBam() writes test.bam alongside a test.bam.bai index
    Rsamtools::asBam(sam, file.path(dir, "test"), overwrite = TRUE)
}

# Both helpers expect whole-genome seqlengths, matching the synthetic header.
testChrLengths <- c(chr1 = 10000, chr2 = 10000)


test_that("fragments import without GenomicRanges attached (#81)", {

    # The #81 crash - MEDIPS raising `could not find function "strand<-"` -
    # only surfaced when GenomicRanges was absent from the search path. The
    # equivalent assertion in test-makeQset.R sits behind skip_long_checks(),
    # which tests/testthat.R enables unconditionally, so it never runs under
    # R CMD check. This block carries no skip_long_checks(), so the guarantee
    # is checked on every run.
    #
    # The precondition is not order-independent: running the suite with
    # skip_long_checks disabled attaches GenomicRanges in an earlier file.
    # Skip rather than report a false failure in that configuration.
    skip_if(
        "package:GenomicRanges" %in% search(),
        "GenomicRanges already attached by an earlier test"
    )
    expect_false("package:GenomicRanges" %in% search())

    dir <- tempfile("mesaBam")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bam <- writeTestBam(dir)

    expect_no_error(
        single <- readSingleEndFragments(
            file = bam, chr.select = "chr1", chr.lengths = testChrLengths
        )
    )
    expect_no_error(
        paired <- readPairedFragments(
            file = bam, chr.select = "chr2", chr.lengths = testChrLengths
        )
    )

    expect_gt(length(single), 0L)
    expect_gt(length(paired), 0L)

    # and the helpers must not have attached it as a side effect
    expect_false("package:GenomicRanges" %in% search())
})


test_that("readSingleEndFragments applies the MEDIPS read filters", {

    dir <- tempfile("mesaBam")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bam <- writeTestBam(dir)

    reads <- readSingleEndFragments(
        file = bam, chr.select = "chr1", chr.lengths = testChrLengths
    )

    # se_long, se_short, se_dupA, se_dupB, se_rev, se_softclip - only the
    # secondary alignment is dropped, and chr2 is excluded.
    expect_equal(length(reads), 6L)
    expect_equal(
        as.character(unique(GenomeInfoDb::seqnames(reads))), "chr1"
    )
    # secondary alignment at 4000 is dropped
    expect_false(4000 %in% BiocGenerics::start(reads))
    # soft-clipped read at 5000 is kept: mesa overrides the MEDIPS
    # simpleCigar default, and has done since f28d678
    expect_true(5000 %in% BiocGenerics::start(reads))
    # strand is stripped only after filtering and deduplication
    expect_true(all(as.character(BiocGenerics::strand(reads)) == "*"))
})


test_that("extend lengthens short reads but never truncates long ones", {

    dir <- tempfile("mesaBam")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bam <- writeTestBam(dir)

    reads <- readSingleEndFragments(
        file = bam, chr.select = "chr1", chr.lengths = testChrLengths,
        extend = 50
    )

    widths <- BiocGenerics::width(reads)
    names(widths) <- BiocGenerics::start(reads)

    # MEDIPS::adjustReads() clamped the extension at pmax(0, extend - width),
    # so the 100 bp read is untouched and the 30 bp read grows to 50.
    expect_equal(unname(widths["1000"]), 100L)
    expect_equal(unname(widths["2000"]), 50L)
    expect_true(all(widths >= 50L))
})


test_that("uniq reproduces the four MEDIPS duplicate-handling branches", {

    dir <- tempfile("mesaBam")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bam <- writeTestBam(dir)

    readChr1 <- function(...) {
        readSingleEndFragments(
            file = bam, chr.select = "chr1",
            chr.lengths = testChrLengths, ...
        )
    }

    expect_equal(length(readChr1(uniq = 0)), 6L)

    # se_dupA and se_dupB collapse; se_rev shares their coordinates but sits
    # on the minus strand, so MEDIPS kept it as a separate location.
    expect_equal(length(readChr1(uniq = 1)), 5L)

    # a p-value caps duplicates per location; at this read depth the Poisson
    # quantile floors to 1, matching uniq = 1
    expect_equal(length(readChr1(uniq = 1e-3)), 5L)

    expect_error(readChr1(uniq = TRUE), "not logical")
    expect_error(readChr1(uniq = 2), "in \\[0, 1\\]")
})


test_that("readPairedFragments reconstructs fragment spans", {

    dir <- tempfile("mesaBam")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bam <- writeTestBam(dir)

    fragments <- readPairedFragments(
        file = bam, chr.select = "chr2", chr.lengths = testChrLengths
    )

    # first mates of pe1, pe_dupA and pe_dupB; the secondary pair is dropped
    expect_equal(length(fragments), 3L)
    expect_false(2000 %in% BiocGenerics::start(fragments))

    # span is the leftmost of read and mate, extended by abs(isize)
    expect_true(all(BiocGenerics::width(fragments) == 200L))
    expect_equal(min(BiocGenerics::start(fragments)), 1000L)
    expect_equal(max(BiocGenerics::end(fragments)), 3199L)

    # the two identical fragments collapse under uniq = 1
    deduped <- readPairedFragments(
        file = bam, chr.select = "chr2", chr.lengths = testChrLengths,
        uniq = 1
    )
    expect_equal(length(deduped), 2L)
})


test_that("chr.select works on a BAM with no index", {

    dir <- tempfile("mesaBam")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bam <- writeTestBam(dir)

    # MEDIPS scanned the whole file and filtered chromosomes afterwards when
    # no index was present; scanBam() would otherwise reject `which`.
    unindexed <- file.path(dir, "noindex.bam")
    file.copy(bam, unindexed)
    expect_false(file.exists(paste0(unindexed, ".bai")))

    expect_message(
        reads <- readSingleEndFragments(
            file = unindexed, chr.select = "chr1",
            chr.lengths = testChrLengths
        ),
        "No BAM index"
    )

    expect_equal(length(reads), 6L)
    expect_equal(
        as.character(unique(GenomeInfoDb::seqnames(reads))), "chr1"
    )
})


test_that("getCGPositions returns one-base motif positions", {

    skip_if_not_installed("BSgenome.Scerevisiae.UCSC.sacCer3")

    positions <- getCGPositions(
        "BSgenome.Scerevisiae.UCSC.sacCer3", chr.select = "chrM"
    )

    # MEDIPS.getPositions() returned IRanges(start, end = start). Width 2
    # would make filter_by_non_overlaps() treat a read starting on the G of a
    # CpG as containing "CG".
    expect_true(all(BiocGenerics::width(positions) == 1L))
    expect_gt(length(positions), 0L)
})
