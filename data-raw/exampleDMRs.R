## code to prepare `exampleDMRs` dataset goes here

## Needs (not in DESCRIPTION): usethis, TxDb.Hsapiens.UCSC.hg38.knownGene,
## org.Hs.eg.db.

## Example annotated set of DMRs generated from the example tumour/normal
## qseaSet, with every pairwise contrast of the `type` column.

data("exampleTumourNormal", package = "mesa")

exampleDMRs <- exampleTumourNormal %>%
  calculateDMRs(variable = "type", contrasts = "all") %>%
  annotateWindows(
    TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    annoDb = "org.Hs.eg.db"
  )

usethis::use_data(exampleDMRs, overwrite = TRUE, compress = "xz")
