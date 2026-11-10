## code to prepare `exampleMouse` dataset goes here

## Needs (not in DESCRIPTION): usethis.

# Example qseaSet data: a small region of chromosome 5 for four mouse
# hematopoietic stem cell (HSC) samples; two replicates of two conditions.
exampleMouse <- readr::read_rds("<PATH TO M073>") %>% 
  filterWindows(seqnames == "chr5", end < 143009500, start > 142500000) %>% 
  filter(stringr::str_detect(sample_name,"Day1|Day5")) %>%
  renameQsetNames("Day1","A") %>%
  renameQsetNames("Day5","B") %>%
  mutate(group = stringr::str_remove(sample_name, "_Rep[12]$")) %>%
  dplyr::select(-project_name) %>%
  # No Mart is stored: exampleMouse is mm10 (GRCm38), and no live Ensembl
  # archive serves that build, so users pass a `mart` to plotGeneHeatmap().
  setMart(mart = NULL)

usethis::use_data(exampleMouse, overwrite = TRUE, compress = "xz")
