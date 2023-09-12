## code to prepare `BSgenome.Mmusculus.UCSC.mm10.CpG.distribution` dataset goes here
BSgenome.Mmusculus.UCSC.mm10.CpG.distribution <- calculateGenomicCGDistribution("BSgenome.Mmusculus.UCSC.mm10")

usethis::use_data(BSgenome.Mmusculus.UCSC.mm10.CpG.distribution, overwrite = TRUE)
