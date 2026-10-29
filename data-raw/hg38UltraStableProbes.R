## code to prepare `ultraStableGR` dataset goes here

## Ultra-stable CpGs: Additional file 2 of Edgar et al. 2014, Epigenetics &
## Chromatin 7:28 (doi:10.1186/1756-8935-7-28), and the InfiniumAnnotation hm450
## hg38 manifest (Zhou, Laird and Shen 2017, doi:10.1093/nar/gkw967), the
## archived 202209 release, because the current one has renamed probeID to
## Probe_ID. Both downloaded on 2026-10-29.
## Originally built from the CRUK internal files
## /data/cep/Methylation/refData/Edgar2014_HyperStable_CpGs.txt and
## /data/cep/Methylation/refData/hm450.hg38.manifest.tsv.
## Re-run 2026-10-29: the same data as the shipped object (identical() is
## FALSE only because of the Seqinfo class attribute, which varies with the
## Bioconductor version).

edgarUrl <- paste0("https://static-content.springer.com/esm/",
    "art%3A10.1186%2F1756-8935-7-28/MediaObjects/",
    "13072_2014_333_MOESM2_ESM.txt")
edgarPath <- tempfile(fileext = ".txt")
download.file(edgarUrl, edgarPath, mode = "wb")

manifestUrl <- paste0("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/",
    "78c0c5d23ca6863c4ab26f5e80f31c0e5d05743f/Anno/HM450/archive/202209/",
    "HM450.hg38.manifest.tsv.gz")
manifestPath <- tempfile(fileext = ".tsv.gz")
download.file(manifestUrl, manifestPath, mode = "wb")

hyperStableCpGs <- read_tsv(edgarPath,
                            col_types = cols(Chromosome_37 = "c"))

# Read in the hg38 information for the 450k microarrays, which includes remapped positions
# Make into a GRanges object with just the 2bp of the CpG island.
hg38MethyProbePositionsStable <- read_tsv(manifestPath) %>%
  mutate(CpG_chrm = str_remove(CpG_chrm,"chr")) %>%
  dplyr::filter(CpG_chrm %in% 1:22) %>%
  dplyr::select(c("CpG_chrm","CpG_beg","CpG_end","probeID")) %>%
  dplyr::rename(chr = CpG_chrm, start= CpG_beg, end = CpG_end, ID = probeID) %>%
  mutate(seqnames = as.factor(chr), start = as.numeric(start), end = as.numeric(end)) %>%
  left_join(hyperStableCpGs, ., by = c("Probe_ID" = "ID"))

hg38UltraStableProbes <- hg38MethyProbePositionsStable %>%
  filter(State == "Ultra-stable Methylated") %>%
  dplyr::select(seqnames, start, end, Probe_ID) %>%
  plyranges::as_granges()


usethis::use_data(hg38UltraStableProbes, overwrite = TRUE)
