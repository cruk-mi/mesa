## Needs (not in DESCRIPTION): usethis.

## Example annotated set of DMRs generated from the example tumour/normal qseaSet

exampleDMRs <- exampleTumourNormal %>% 
  calculateDMRs(variable = "type", contrasts = "all") %>% 
  annotateWindows(
    TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene", 
    annoDb = "org.Hs.eg.db"
  )

usethis::use_data(exampleDMRs, overwrite = TRUE, compress = "xz")
