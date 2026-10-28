# Synthetic BAMs for the makeQset() option tests, so they need neither the
# MEDIPSData reads nor a human BSgenome. The reads sit on sacCer3's chrI, and
# each file takes about 0.05 s to write.

# Writes a coordinate-sorted, indexed BAM of `nPairs` proper pairs (50 bp
# reads, MAPQ 60, inserts of 150-300 bp) spread evenly along chrI. Positions
# come from modular arithmetic, not the RNG, so the files never move the
# test's random seed. Returns the BAM path.
writeYeastPairs <- function(
    dir, name, nPairs = 5000L, chrLen = 230218L, offset = 0L
) {

    idx <- seq_len(nPairs)
    insert <- 150L + (idx * 31L) %% 151L
    start1 <- 1L + (idx * 7919L + offset) %% (chrLen - 400L)
    start2 <- start1 + insert - 50L
    readName <- paste0("pair", idx)

    rec <- function(flag, pos, pnext, tlen) {
        paste(readName, flag, "chrI", pos, 60, "50M", "=", pnext, tlen,
            strrep("A", 50), strrep("I", 50), sep = "\t")
    }

    recs <- c(
        rec(99, start1, start2, insert),
        rec(147, start2, start1, -insert)
    )
    recs <- recs[order(c(start1, start2))]

    sam <- file.path(dir, paste0(name, ".sam"))
    writeLines(
        c("@HD\tVN:1.6\tSO:coordinate",
            paste0("@SQ\tSN:chrI\tLN:", chrLen),
            recs),
        sam
    )
    Rsamtools::asBam(sam, file.path(dir, name), overwrite = TRUE)
}

# The sampleTable for two samples with their two input libraries (which CNV
# calling and qseaPaired normalisation need), built once per test run in
# tempdir(), which R removes when the session ends.
.mesaTestBams <- new.env(parent = emptyenv())

yeastSampleTable <- function() {

    if (is.null(.mesaTestBams$sampleTable)) {
        dir <- tempfile("mesaYeastPairs")
        dir.create(dir)
        bams <- vapply(seq_len(4), function(i) {
            writeYeastPairs(dir, paste0("sample", i), offset = i * 101L)
        }, character(1))

        .mesaTestBams$sampleTable <- data.frame(
            sample_name = c("Normal1", "Tumour1"),
            group = c("Normal", "Tumour"),
            file_name = unname(bams[1:2]),
            input_file = unname(bams[3:4])
        )
    }
    .mesaTestBams$sampleTable
}
