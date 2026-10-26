## code to prepare `gc_hg19_chr22_100kb` dataset goes here
##
## GC fraction of 100 kb bins on chr22 of hg19, computed from the public
## BSgenome.Hsapiens.UCSC.hg19 sequence. Bins that are all N get -1, the value
## HMMcopy uses for a bin without a GC estimate.

binSize <- 100000L
genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
chrLength <- GenomeInfoDb::seqlengths(genome)[["chr22"]]

# the same tiling as addHMMcopyCNV(), which drops the last partial bin
starts <- seq(1L, chrLength - binSize + 1L, by = binSize)
bins <- Biostrings::Views(genome[["chr22"]], starts, starts + binSize - 1L)
counts <- Biostrings::letterFrequency(bins, c("G", "C", "A", "T"))

gc_hg19_chr22_100kb <- data.frame(
    chr = "chr22",
    start = starts,
    end = starts + binSize - 1L,
    gc = rowSums(counts[, c("G", "C")]) / rowSums(counts)
)
gc_hg19_chr22_100kb$gc[is.na(gc_hg19_chr22_100kb$gc)] <- -1

usethis::use_data(gc_hg19_chr22_100kb, overwrite = TRUE, compress = "xz")
