library(tidyverse)

lines_df <- readLines("DESCRIPTION") %>%
  enframe()

import_line <- lines_df %>%
  filter(value == "Imports:") %>%
  pull(name)

biocViews_line <- lines_df %>%
  filter(str_detect(value,"biocViews:")) %>%
  pull(name)

requiredPackages <- lines_df %>%
  filter(name > import_line, name < biocViews_line) %>%
  filter(!str_detect(value,"Depends:|^R\\(")) %>%
  mutate(value = str_remove_all(value, " |,"),
         value = str_remove_all(value,"\\(.*|VignetteBuilder:")) %>%
  filter(value != "R") %>%
  pull(value)

remotes::install_cran(requiredPackages, repos = BiocManager::repositories())