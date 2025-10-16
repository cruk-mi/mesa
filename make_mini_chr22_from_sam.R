#!/usr/bin/env Rscript
# Create inst/extdata/mini_chr22.{bam,bai} from a minimal SAM on chr22

message("→ Picking installed human BSgenome…")

# pick whichever genome package is installed
genome_pkg <- if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
    "BSgenome.Hsapiens.UCSC.hg19"
} else if (requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
    "BSgenome.Hsapiens.NCBI.GRCh38"
} else {
    stop("Install BSgenome.Hsapiens.UCSC.hg19 or BSgenome.Hsapiens.NCBI.GRCh38")
}

if (!requireNamespace("GenomeInfoDb", quietly = TRUE))
    stop("Please install GenomeInfoDb")

# load BSgenome object
genome <- getExportedValue(genome_pkg, genome_pkg)

# identify build + seqstyle (no custom %||%)
genome_build <- NA
try({ genome_build <- GenomeInfoDb::genome(genome)[1] }, silent = TRUE)
if (is.null(genome_build) || is.na(genome_build) || length(genome_build) == 0) genome_build <- "unknown"
seqstyle <- if ("chr1" %in% GenomeInfoDb::seqlevels(genome)) "UCSC-style (chr1..chr22)" else "NCBI-style (1..22)"
message(sprintf("✔ Using %s [%s, %s]", genome_pkg, genome_build, seqstyle))

# choose chr22 (or 22)
chr <- if ("chr22" %in% GenomeInfoDb::seqlevels(genome)) "chr22" else "22"
len <- GenomeInfoDb::seqlengths(genome)[chr]
if (is.na(len)) stop("Couldn't get length for ", chr)

# output locations
outdir <- file.path("inst", "extdata")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sam_path  <- file.path(outdir, "mini_chr22.sam")
dest_base <- sub("\\.sam$", "", sam_path)   # asBam destination *without* .bam
bam_path  <- paste0(dest_base, ".bam")

# design: one proper pair (99/147) + one unpaired R1 (0)
start0 <- 20000000L
readL  <- 100L
gap    <- 200L
pos1   <- start0 + 100L
pos2   <- pos1 + gap
seq100 <- paste(rep("A", readL), collapse = "")
qual100<- paste(rep("I", readL), collapse = "")

hdr <- sprintf("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:%s\tLN:%d", chr, as.integer(len))
sam_lines <- c(
    hdr,
    sprintf("r001\t99\t%s\t%d\t60\t%dM\t=\t%d\t%d\t%s\t%s",
            chr, pos1, readL, pos2,  gap + readL, seq100, qual100),
    sprintf("r001\t147\t%s\t%d\t60\t%dM\t=\t%d\t-%d\t%s\t%s",
            chr, pos2, readL, pos1,  gap + readL, seq100, qual100),
    sprintf("r002\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s",
            chr, start0 + 500L, readL, seq100, qual100)
)

writeLines(sam_lines, sam_path)
message("✓ Wrote ", sam_path)

if (!requireNamespace("Rsamtools", quietly = TRUE))
    stop("Please install Rsamtools")

Rsamtools::asBam(sam_path, destination = dest_base, overwrite = TRUE)
Rsamtools::indexBam(bam_path)
unlink(sam_path)

message("✅ Created: ", bam_path, " and ", bam_path, ".bai")
