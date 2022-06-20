## code to prepare `ultraStableGR` dataset goes here

hyperStableCpGs <- read_tsv("/data/cep/Methylation/refData/Edgar2014_HyperStable_CpGs.txt",
                            col_types = cols(Chromosome_37 = "c"))

# Read in the hg38 information for the 450k microarrays, which includes remapped positions
# Make into a GRanges object with just the 2bp of the CpG island.
hg38MethyProbePositionsStable <- read_tsv("/data/cep/Methylation/refData/hm450.hg38.manifest.tsv") %>%
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
