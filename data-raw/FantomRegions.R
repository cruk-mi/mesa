## code to prepare `FantomRegions` dataset goes here

FantomRegions <- plyranges::read_bed("/data/cep/Methylation/refData/FANTOM5_hg38_noChr.bed") %>%
  as_tibble() %>%
  dplyr::select(seqnames, start, end)

usethis::use_data(FantomRegions, overwrite = TRUE)
