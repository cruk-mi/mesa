## code to prepare `map_hg38_10kb` dataset goes here

map_hg38_10kb <- HMMcopy::wigToRangedData("/data/cep/Methylation/refData/map_hg38_10kb.wig", verbose = FALSE) %>%
  dplyr::mutate(chr = stringr::str_remove(chr,"chr"), end = end - 1) %>%
  dplyr::rename(map = value)

usethis::use_data(map_hg38_10kb, overwrite = TRUE, compress = "xz", compression_level = 9)
