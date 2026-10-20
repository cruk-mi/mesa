## code to prepare `exampleMouse` dataset goes here

# M073 mouse samples
exampleMouse <- read_rds("<PATH TO M073>") %>% 
  filterWindows(seqnames == "chr5", end < 143009500, start > 142500000) %>% 
  filter(str_detect(sample_name,"Day1|Day5")) %>%
  renameQsetNames("Day1","A") %>%
  renameQsetNames("Day5","B") %>%
  # No Mart is stored: exampleMouse is mm10 (GRCm38), and no live Ensembl
  # archive serves that build, so users pass a `mart` to plotGeneHeatmap().
  setMart(mart = NULL)

usethis::use_data(exampleMouse, overwrite = TRUE, compress = "xz")
