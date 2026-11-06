## code to prepare `DATASET` dataset goes here

## Needs (not in DESCRIPTION): usethis.

BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution <- calculateGenomicCGDistribution("BSgenome.Hsapiens.NCBI.GRCh38")

usethis::use_data(BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution, overwrite = TRUE)
