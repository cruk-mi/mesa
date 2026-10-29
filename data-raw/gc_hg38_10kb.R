## code to prepare `gc_hg38_10kb` dataset goes here

## GC content of 10 kb bins of GRCh38 (HMMcopy wig format), from the
## ichorCNA repository (https://github.com/broadinstitute/ichorCNA), the source
## named in the help page. Downloaded on 2026-10-29.
## Originally built from the CRUK internal file
## /data/cep/Methylation/refData/gc_hg38_10kb.wig.
## Re-run 2026-10-29: no shipped object to compare with.

url <- paste0("https://github.com/broadinstitute/ichorCNA/raw/",
    "5bfc03ed854f0e93fe5b624c97c1290fa0053837/inst/extdata/",
    "gc_hg38_10kb.wig")
path <- tempfile(fileext = ".wig")
download.file(url, path, mode = "wb")

gc_hg38_10kb <- HMMcopy::wigToRangedData(path, verbose = FALSE) %>%
  dplyr::mutate(chr = stringr::str_remove(chr,"chr"), end = end - 1) %>%
  dplyr::rename(gc = value)

usethis::use_data(gc_hg38_10kb, overwrite = TRUE, compress = "xz")
