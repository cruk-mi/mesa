## code to prepare `map_hg38_1000kb` dataset goes here

map_hg38_50kb <- HMMcopy::wigToRangedData("/data/cep/Methylation/refData/map_hg38_50kb.wig", verbose = FALSE) %>%
  dplyr::mutate(chr = stringr::str_remove(chr,"chr"), end = end - 1) %>%
  dplyr::rename(map = value)

usethis::use_data(map_hg38_50kb, overwrite = TRUE)
