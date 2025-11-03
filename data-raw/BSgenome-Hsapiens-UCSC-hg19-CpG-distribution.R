## code to prepare `BSgenome.Hsapiens.UCSC.hg19.CpG.distribution` dataset goes here

BSgenome.Hsapiens.UCSC.hg19.CpG.distribution <- calculateGenomicCGDistribution("BSgenome.Hsapiens.UCSC.hg19")

usethis::use_data(BSgenome.Hsapiens.UCSC.hg19.CpG.distribution, overwrite = TRUE)
