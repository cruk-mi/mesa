remotes::install_cran(c("tibble","dplyr","stringr"))

lines_df <- readLines("DESCRIPTION") |>
  tibble::enframe()

import_line <- lines_df |>
  dplyr::filter(value == "Imports:") |>
  dplyr::pull(name)

biocViews_line <- lines_df |>
  dplyr::filter(stringr::str_detect(value,"biocViews:")) |>
  dplyr::pull(name)

requiredPackages <- lines_df |>
  dplyr::filter(name > import_line, name < biocViews_line) |>
  dplyr::filter(!stringr::str_detect(value,"Depends:|^R\\(")) |>
  dplyr::mutate(value = stringr::str_remove_all(value, " |,"),
         value = stringr::str_remove_all(value,"\\(.*|VignetteBuilder:")) |>
  dplyr::filter(value != "R") |>
  dplyr::pull(value)

remotes::install_cran(requiredPackages, repos = BiocManager::repositories())