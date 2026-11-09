## code to prepare `hg38UltraStableProbes` dataset goes here

## Needs (not in DESCRIPTION): usethis.

## Ultra-stable CpGs: Additional file 2 of Edgar et al. 2014, Epigenetics &
## Chromatin 7:28 (doi:10.1186/1756-8935-7-28), downloaded on 2026-10-29.
##
## Depends on the package data object `hg38_450kArrayGR` (made by
## hg38_450kArrayData.R) for the hg38 probe positions. Re-run that script
## first if the object changes.

edgarUrl <- paste0("https://static-content.springer.com/esm/",
    "art%3A10.1186%2F1756-8935-7-28/MediaObjects/",
    "13072_2014_333_MOESM2_ESM.txt")
edgarPath <- tempfile(fileext = ".txt")
download.file(edgarUrl, edgarPath, mode = "wb")

hyperStableCpGs <- readr::read_tsv(edgarPath,
                                   col_types = readr::cols(Chromosome_37 = "c"))

# Positions of the 450k probes on the autosomes, from the package data object.
hg38MethyProbePositionsStable <- as.data.frame(mesa::hg38_450kArrayGR) %>%
  dplyr::filter(seqnames %in% 1:22) %>%
  dplyr::transmute(seqnames = as.factor(as.character(seqnames)),
                   start, end, ID) %>%
  dplyr::left_join(hyperStableCpGs, ., by = c("Probe_ID" = "ID"))

hg38UltraStableProbes <- hg38MethyProbePositionsStable %>%
  dplyr::filter(State == "Ultra-stable Methylated") %>%
  dplyr::select(seqnames, start, end, Probe_ID) %>%
  plyranges::as_granges()

usethis::use_data(hg38UltraStableProbes, overwrite = TRUE)
