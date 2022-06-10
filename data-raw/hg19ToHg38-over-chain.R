## code to prepare `hg19ToHg38.over.chain` dataset goes here

#downloaded from ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

hg19ToHg38.over.chain <- rtracklayer::import.chain("/data/cep/Methylation/refData/hg19ToHg38.over.chain")

usethis::use_data(hg19ToHg38.over.chain, overwrite = TRUE)
