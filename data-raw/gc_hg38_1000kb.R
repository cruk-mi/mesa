## code to prepare `gc_hg38_1000kb.wig` dataset goes here

## GC content of 1 Mb bins of GRCh38 (HMMcopy wig format), from the
## ichorCNA repository (https://github.com/broadinstitute/ichorCNA), the source
## named in the help page. Downloaded on 2026-10-29.

url <- paste0("https://github.com/broadinstitute/ichorCNA/raw/master/",
    "inst/extdata/gc_hg38_1000kb.wig")
path <- tempfile(fileext = ".wig")
download.file(url, path, mode = "wb")

gc_hg38_1000kb <- HMMcopy::wigToRangedData(path, verbose = FALSE) %>%
  dplyr::mutate(chr = stringr::str_remove(chr,"chr"), end = end - 1) %>%
  dplyr::rename(gc = value)

usethis::use_data(gc_hg38_1000kb, overwrite = TRUE, compress = "xz")
