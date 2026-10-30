## code to prepare `gc_hg38_1000kb.wig` dataset goes here

## GC content of 1 Mb bins of GRCh38 (HMMcopy wig format), from the
## ichorCNA repository (https://github.com/broadinstitute/ichorCNA).

url <- paste0("https://github.com/broadinstitute/ichorCNA/raw/",
    "5bfc03ed854f0e93fe5b624c97c1290fa0053837/inst/extdata/",
    "gc_hg38_1000kb.wig")
path <- tempfile(fileext = ".wig")
download.file(url, path, mode = "wb")

gc_hg38_1000kb <- HMMcopy::wigToRangedData(path, verbose = FALSE) %>%
  dplyr::mutate(chr = stringr::str_remove(chr,"chr"), end = end - 1) %>%
  dplyr::rename(gc = value)

usethis::use_data(gc_hg38_1000kb, overwrite = TRUE, compress = "xz")
