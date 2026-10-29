## code to prepare `map_hg38_1000kb` dataset goes here

## Mappability of 1 Mb bins of GRCh38 (HMMcopy wig format), from the
## ichorCNA repository (https://github.com/broadinstitute/ichorCNA), the source
## named in the help page. Downloaded on 2026-10-29.

url <- paste0("https://github.com/broadinstitute/ichorCNA/raw/master/",
    "inst/extdata/map_hg38_1000kb.wig")
path <- tempfile(fileext = ".wig")
download.file(url, path, mode = "wb")

map_hg38_1000kb <- HMMcopy::wigToRangedData(path, verbose = FALSE) %>%
  dplyr::mutate(chr = stringr::str_remove(chr,"chr"), end = end - 1) %>%
  dplyr::rename(map = value)

usethis::use_data(map_hg38_1000kb, overwrite = TRUE, compress = "xz")
