## code to prepare `hg19ToHg38.over.chain` dataset goes here

## UCSC liftOver chain from hg19 to hg38, downloaded on 2026-10-29. The file
## name hg19ToHg38 is the direction the object lifts (hg19 -> hg38); the old
## comment gave the opposite hg38ToHg19 file.
## Originally built from the CRUK internal file
## /data/cep/Methylation/refData/hg19ToHg38.over.chain.
## Re-run 2026-10-29: the same as the shipped object.

url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/",
    "hg19ToHg38.over.chain.gz")
path <- tempfile(fileext = ".over.chain.gz")
download.file(url, path, mode = "wb")

chainFile <- tempfile(fileext = ".over.chain")
writeLines(readLines(path), chainFile)

hg19ToHg38.over.chain <- rtracklayer::import.chain(chainFile)

usethis::use_data(hg19ToHg38.over.chain, overwrite = TRUE, compress = "xz")
