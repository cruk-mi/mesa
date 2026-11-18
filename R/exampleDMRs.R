#' Annotated example differentially methylated regions
#'
#' DMR results for the tumour-versus-normal comparison in
#' [exampleTumourNormal], annotated with hg38 gene information.
#'
#' @format A tibble with 23 rows and 24 columns, containing DMR statistics and
#' gene annotations.
#'
#' @source Generated with [calculateDMRs()] and [annotateWindows()] using
#' `TxDb.Hsapiens.UCSC.hg38.knownGene` and `org.Hs.eg.db`.
#'
#' @examples
#' data("exampleDMRs", package = "mesa")
#' head(exampleDMRs)
"exampleDMRs"
