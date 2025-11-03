## code to prepare `fantomPlusRegions` dataset goes here

fantomPlusRegions <- hg38CpGIslands %>%
       plyranges::anchor_centre() %>%
       plyranges::stretch(8000) %>%
       plyranges::union_ranges(FantomRegions) %>%
       plyranges::reduce_ranges()

usethis::use_data(fantomPlusRegions, overwrite = TRUE)
