#!/usr/bin/env Rscript
# Create inst/extdata/mini_chr22_GRCh38.{bam,bai} using NCBI GRCh38 contig style (no 'chr' prefix)

message("→ Using NCBI GRCh38 (no 'chr' prefix)")

# Require GRCh38 BSgenome
if (!requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
    stop("Please install BSgenome.Hsapiens.NCBI.GRCh38")
}
if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop("Please install GenomeInfoDb")
}

# Load BSgenome object
genome_pkg <- "BSgenome.Hsapiens.NCBI.GRCh38"
genome <- getExportedValue(genome_pkg, genome_pkg)

# Sanity + info
if (!"22" %in% GenomeInfoDb::seqlevels(genome)) {
    stop("GRCh38 genome detected but seqlevel '22' not found. Check your BSgenome installation.")
}
genome_build <- tryCatch(GenomeInfoDb::genome(genome)[1], error = function(e) NA)
if (is.null(genome_build) || is.na(genome_build) || length(genome_build) == 0) genome_build <- "GRCh38"
message(sprintf("✔ Using %s [%s, NCBI-style (1..22)]", genome_pkg, genome_build))

# Target contig and length (NCBI style)
chr <- "22"
len <- GenomeInfoDb::seqlengths(genome)[chr]
if (is.na(len)) stop("Couldn't get length for ", chr)

# Output paths
outdir <- file.path("inst", "extdata")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sam_path  <- file.path(outdir, "mini_chr22_GRCh38.sam")   # keep filename for compatibility
dest_base <- sub("\\.sam$", "", sam_path)          # asBam destination (no .bam)
bam_path  <- paste0(dest_base, ".bam")

# Design: one proper pair (flags 99/147) + one unpaired R1 (flag 0)
start0 <- 20000000L
readL  <- 100L
gap    <- 200L
pos1   <- start0 + 100L
pos2   <- pos1 + gap
seq100 <- paste(rep("A", readL), collapse = "")
qual100<- paste(rep("I", readL), collapse = "")

# SAM header: SN:22 (NCBI style)
hdr <- sprintf("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:%s\tLN:%d", chr, as.integer(len))
sam_lines <- c(
    hdr,
    # qname  flag RNAME POS    MAPQ CIGAR RNEXT PNEXT TLEN     SEQ    QUAL
    sprintf("r001\t99\t%s\t%d\t60\t%dM\t=\t%d\t%d\t%s\t%s",
            chr, pos1, readL, pos2,  gap + readL, seq100, qual100),
    sprintf("r001\t147\t%s\t%d\t60\t%dM\t=\t%d\t-%d\t%s\t%s",
            chr, pos2, readL, pos1,  gap + readL, seq100, qual100),
    sprintf("r002\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s",
            chr, start0 + 500L, readL,           seq100, qual100)
)

writeLines(sam_lines, sam_path)
message("✓ Wrote ", sam_path)

# Convert to BAM + index
if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Please install Rsamtools")
}
Rsamtools::asBam(sam_path, destination = dest_base, overwrite = TRUE)
Rsamtools::indexBam(bam_path)
unlink(sam_path)

message("✅ Created: ", bam_path, " and ", bam_path, ".bai")
