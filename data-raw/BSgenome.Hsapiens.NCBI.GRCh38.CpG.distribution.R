## code to prepare `DATASET` dataset goes here

BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution <- calculateGenomicCGDistribution("BSgenome.Hsapiens.NCBI.GRCh38")

usethis::use_data(BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution, overwrite = TRUE)
