#' Annotated example differentially methylated regions
#'
#' DMR results for all pairwise contrasts between the `type` groups of
#' [exampleTumourNormal], annotated with hg38 gene information.
#'
#' @format A tibble with 109 rows and 57 columns, in the wide format returned
#' by [calculateDMRs()]: one set of statistic columns per contrast (ten
#' contrasts between NormalColon, CRC, NormalLung, LUAD and LUSC), the
#' per-group `_beta_means` and `_nrpm_means` columns, and the gene annotation
#' columns added by [annotateWindows()]. Use [pivotDMRsLonger()] to convert
#' it to long format.
#'
#' @source Generated from [exampleTumourNormal] with [calculateDMRs()] (all
#' contrasts of `type`) and [annotateWindows()], using
#' `TxDb.Hsapiens.UCSC.hg38.knownGene` and `org.Hs.eg.db`. See
#' `data-raw/exampleDMRs.R`.
#'
#' @examples
#' data("exampleDMRs", package = "mesa")
#' head(exampleDMRs)
"exampleDMRs"
