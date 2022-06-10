## code to prepare `GRCh38.99.gtf` dataset goes here

GRCh38.99.gtf <- rtracklayer::import("/mnt/gpfs2/data/cep/CDX_Model_NGS_Data/Reference/Homo_sapiens.GRCh38.99.gtf") %>%
  filter(type %in% c("gene", "exon","five_prime_utr","three_prime_utr")) %>%
  dplyr::select(type, gene_id,gene_name)

usethis::use_data(GRCh38.99.gtf, overwrite = TRUE)
