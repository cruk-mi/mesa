## code to prepare `FantomRegions` dataset goes here

## FANTOM5 hg38 enhancers (Andersson et al. 2014, doi:10.1038/nature12787),
## reprocessed to hg38 by FANTOM, downloaded on 2026-10-29. Remove "chr".
## Originally built from the CRUK internal file
## /data/cep/Methylation/refData/FANTOM5_hg38_noChr.bed.
## Re-run 2026-10-29: the same as the shipped object.

url <- paste0("https://fantom.gsc.riken.jp/5/datafiles/reprocessed/",
    "hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz")
path <- tempfile(fileext = ".bed.gz")
download.file(url, path, mode = "wb")

bedFile <- tempfile(fileext = ".bed")
writeLines(sub("^chr", "", readLines(path)), bedFile)

FantomRegions <- plyranges::read_bed(bedFile) %>%
  as_tibble() %>%
  dplyr::select(seqnames, start, end)

usethis::use_data(FantomRegions, overwrite = TRUE)
