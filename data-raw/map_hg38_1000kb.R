## code to prepare `map_hg38_1000kb` dataset goes here

map_hg38_1000kb <- HMMcopy::wigToRangedData("/data/cep/Methylation/refData/map_hg38_1000kb.wig", verbose = FALSE) %>%
  dplyr::mutate(chr = stringr::str_remove(chr,"chr"), end = end - 1) %>%
  dplyr::rename(map = value)

usethis::use_data(map_hg38_1000kb, overwrite = TRUE, compress = "xz")
