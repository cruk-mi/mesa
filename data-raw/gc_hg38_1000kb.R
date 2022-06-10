## code to prepare `gc_hg38_1000kb.wig` dataset goes here

gc_hg38_1000kb <- HMMcopy::wigToRangedData("/data/cep/Methylation/refData/gc_hg38_1000kb.wig", verbose = FALSE) %>%
  dplyr::mutate(chr = stringr::str_remove(chr,"chr"), end = end - 1) %>%
  dplyr::rename(gc = value)

usethis::use_data(gc_hg38_1000kb, overwrite = TRUE)
