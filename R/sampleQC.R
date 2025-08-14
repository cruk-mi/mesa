#' Summarise the fraction of hyper-stable methylated regions (GRCh38 only)
#'
#' Calculates, for each sample in a `qseaSet`, the proportion of *Edgar et al.* “ultra-stable”
#' regions that are present and clearly methylated (beta ≥ `minBeta`) after filtering by CpG
#' density (≥ `minDensity`). The result is added as the column `hyperStableEdgar` in the
#' sample table. If that column already exists it will be overwritten.
#'
#' @details
#' - **Genome restriction:** Only defined for `BSgenome.Hsapiens.NCBI.GRCh38`; the function stops otherwise.
#' - **Regions used:** Requires `mesa::hg38UltraStableProbes` (GRCh38 “ultra-stable” methylated windows
#'   derived from array data).
#' - **Computation:** Overlap windows with `hg38UltraStableProbes`, build a beta table
#'   (`qsea::makeTable(norm_methods = "beta")`), keep windows with `CpG_density >= minDensity`,
#'   treat `NA` beta values as 0 for thresholding, and take the fraction ≥ `minBeta` per sample.
#'
#' @param qseaSet A `qseaSet` object (must be GRCh38).
#' @param minDensity Numeric. Minimum `CpG_density` to keep a window (default `5`).
#' @param minBeta Numeric. Beta threshold to call a window methylated (default `0.8`).
#'
#' @return
#' The input `qseaSet`, with its sample table augmented by a numeric column `hyperStableEdgar`
#' (range 0–1) giving, per sample, the fraction of hyper-stable regions called methylated.
#' The `qseaSet` is returned.
#'
#' @seealso
#' [qsea::makeTable()], [qsea::getSampleTable()], [qsea::getSampleNames()], [hg38UltraStableProbes],
#' [getSampleQCSummary()]
#'
#' @references
#' Edgar R, Tan PPC, Portales-Casamar E, Pavlidis P (2014).
#' *Meta-analysis of human methylomes reveals stably methylated sequences surrounding CpG islands associated with high gene expression*.
#' <https://pubmed.ncbi.nlm.nih.gov/25493099/>
#'
#' @examples
#' \donttest{
#' if (requireNamespace("qsea", quietly = TRUE)) {
#'   # Example qseaSet shipped with mesa (should be GRCh38)
#'   if (system.file("data", "exampleTumourNormal.rda", package = "mesa") != "") {
#'     data("exampleTumourNormal", package = "mesa")
#'
#'     # Add the hyper-stable fraction with defaults
#'     qset1 <- addHyperStableFraction(exampleTumourNormal)
#'     head(qsea::getSampleTable(qset1)$hyperStableEdgar)
#'
#'     # Stricter thresholds
#'     qset2 <- addHyperStableFraction(exampleTumourNormal, minDensity = 10, minBeta = 0.9)
#'     summary(qsea::getSampleTable(qset2)$hyperStableEdgar)
#'   }
#' }
#' }
#' @export
addHyperStableFraction <- function(qseaSet, minDensity = 5, minBeta = 0.8){

  if(qsea:::getGenome(qseaSet) != "BSgenome.Hsapiens.NCBI.GRCh38") {
    stop("This function is only currently defined for BSgenome.Hsapiens.NCBI.GRCh38.")
  }
  
  if ("hyperStableEdgar" %in% colnames(qsea::getSampleTable(qseaSet))) {
    qseaSet <- qseaSet %>%
      selectQset(-tidyselect::matches("^hyperStableEdgar$"))
  }

  hyperStableBetaTable <- qseaSet %>%
    filterByOverlaps(mesa::hg38UltraStableProbes) %>%
    qsea::makeTable(norm_methods = "beta", samples = qsea::getSampleNames(.)) %>%
    dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
    dplyr::filter(CpG_density >= 5) %>%
    dplyr::select(tidyselect::matches("beta"))

  overp8Fraction <- hyperStableBetaTable %>%
    replace(., is.na(.), 0) %>%
    {. >= 0.8} %>%
    apply(2,mean,na.rm = TRUE)

  newData <- tibble::tibble(sample_name = names(overp8Fraction),
                            hyperStableEdgar = overp8Fraction) %>%
    dplyr::mutate(sample_name = stringr::str_remove(sample_name,"_beta")) %>%
    dplyr::mutate(hyperStableEdgar = tidyr::replace_na(round(hyperStableEdgar, 3),0))

  qseaSet@sampleTable <- cbind(qseaSet@sampleTable, dplyr::select(newData,-sample_name))

  return(qseaSet)

}

#' Summarise key QC fields per sample
#'
#' Returns a tidy summary of the most relevant QC fields from a `qseaSet` sample table.
#' If library metrics are missing, they are added via `addLibraryInformation()` before summarising.
#'
#' @param qseaSet A `qseaSet` object.
#'
#' @return
#' A tibble with one row per sample and (at least) the columns:
#' - `sample_name`
#' - any columns matching `valid_fragment`, `relH`, `hyperStable`, or `ichorTumo` (if present).
#' Rows are ordered by `sample_name`.
#'
#' @seealso
#' [qsea::getSampleTable()], [addLibraryInformation()], [addHyperStableFraction()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("qsea", quietly = TRUE)) {
#'   if (system.file("data", "exampleTumourNormal.rda", package = "mesa") != "") {
#'     data("exampleTumourNormal", package = "mesa")
#'     # Ensure the hyper-stable column exists to demonstrate inclusion
#'     qset <- addHyperStableFraction(exampleTumourNormal)
#'     qc   <- getSampleQCSummary(qset)
#'     qc
#'   }
#' }
#' }
#' @export
getSampleQCSummary <- function(qseaSet){

  if (!("valid_fragments" %in% colnames(qsea::getSampleTable(qseaSet)))) {
    qseaSet <- qseaSet %>%
      addLibraryInformation()
  }

  qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::select(sample_name, tidyselect::matches("valid_fragment|relH|hyperStable|ichorTumo")) %>%
    dplyr::arrange(sample_name) %>%
    tibble::as_tibble() %>%
    return()

}

