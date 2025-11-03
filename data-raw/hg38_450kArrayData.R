## code to prepare `hg38_450kArrayGR` dataset goes here

hg38_450kArrayGR <- readr::read_tsv("/data/cep/Methylation/refData/hm450.hg38.manifest.tsv") %>%
  dplyr::mutate(CpG_chrm = stringr::str_remove(CpG_chrm,"chr")) %>%
  dplyr::filter(CpG_chrm %in% c(1:22,"X","Y")) %>%
  dplyr::select(c("CpG_chrm","CpG_beg","CpG_end","probeID","MASK_general")) %>%
  dplyr::rename(chr = CpG_chrm, start= CpG_beg, end = CpG_end, ID = probeID) %>%
  dplyr::mutate(seqnames = as.factor(chr), start = as.numeric(start), end = as.numeric(end)) %>%
  dplyr::select(-chr) %>%
  plyranges::as_granges()

usethis::use_data(hg38_450kArrayGR, overwrite = TRUE, compress = "xz", compression_level = 9)
