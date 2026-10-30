## code to prepare `hg38_450kArrayGR` dataset goes here

## Illumina HM450 probes remapped to hg38 (Zhou, Laird and Shen 2017,
## doi:10.1093/nar/gkw967), the InfiniumAnnotation hm450 hg38 manifest. This is
## the archived 202209 release: the current manifest has renamed probeID to
## Probe_ID and dropped the MASK_general column.

url <- paste0("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/",
    "78c0c5d23ca6863c4ab26f5e80f31c0e5d05743f/Anno/HM450/archive/202209/",
    "HM450.hg38.manifest.tsv.gz")
manifest <- tempfile(fileext = ".tsv.gz")
download.file(url, manifest, mode = "wb")

hg38_450kArrayGR <- readr::read_tsv(manifest) %>%
  dplyr::mutate(CpG_chrm = stringr::str_remove(CpG_chrm,"chr")) %>%
  dplyr::filter(CpG_chrm %in% c(1:22,"X","Y")) %>%
  dplyr::select(c("CpG_chrm","CpG_beg","CpG_end","probeID","MASK_general")) %>%
  dplyr::rename(chr = CpG_chrm, start= CpG_beg, end = CpG_end, ID = probeID) %>%
  dplyr::mutate(seqnames = as.factor(chr), start = as.numeric(start), end = as.numeric(end)) %>%
  dplyr::select(-chr) %>%
  plyranges::as_granges()

usethis::use_data(hg38_450kArrayGR, overwrite = TRUE, compress = "xz", compression_level = 9)
