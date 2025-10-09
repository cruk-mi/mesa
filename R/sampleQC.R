#' Summarise the fraction of hyper-stable methylated regions (GRCh38 only)
#'
#' For each sample, compute the proportion of *Edgar et al.* “ultra-stable”
#' windows that are clearly methylated (beta ≥ `minBeta`) after filtering by
#' CpG density (≥ `minDensity`). The result is added to the sample table as
#' `hyperStableEdgar` (overwriting any existing column of that name).
#'
#' @details
#' - **Genome restriction.** Only supported for
#'   `"BSgenome.Hsapiens.NCBI.GRCh38"`; the function stops otherwise (checked via
#'   `qsea:::getGenome()`).
#' - **Regions used.** Uses `mesa::hg38UltraStableProbes`, a GRCh38 set derived
#'   from array data (*Edgar et al.*, 2014).
#' - **Computation.** Windows overlapping the ultra-stable set are extracted,
#'   beta values are computed via `qsea::makeTable(norm_methods = "beta")`,
#'   windows with `CpG_density < minDensity` are dropped, `NA` beta are treated
#'   as `0` for thresholding, and the per-sample fraction with beta ≥ `minBeta`
#'   is returned.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object (must be GRCh38).  
#'   **Default:** none (must be supplied).
#'
#' @param minDensity `numeric(1)`.  
#'   Minimum `CpG_density` to keep a window.  
#'   **Default:** `5`.
#'
#' @param minBeta `numeric(1)`.  
#'   Beta threshold to call a window methylated.  
#'   **Default:** `0.8`.
#'
#' @return
#' The input `qseaSet` with its `sampleTable` augmented by a numeric column
#' `hyperStableEdgar` (range 0–1) giving, per sample, the fraction of ultra-stable
#' windows called methylated.
#'
#' @seealso
#' [qsea::makeTable()], [qsea::getSampleTable()], [qsea::getSampleNames()],
#' [hg38UltraStableProbes], [getSampleQCSummary()]
#'
#' @references
#' Edgar R, Tan PPC, Portales-Casamar E, Pavlidis P (2014).  
#' *Meta-analysis of human methylomes reveals stably methylated sequences
#' surrounding CpG islands associated with high gene expression*.  
#' <https://pubmed.ncbi.nlm.nih.gov/25493099/>
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Run only if the example windows overlap the ultra-stable probes
#' exampleTumourNormal %>%
#' {
#'   has_ov <- sum(
#'     GenomicRanges::countOverlaps(
#'       qsea::getRegions(.),
#'       mesa::hg38UltraStableProbes
#'     )
#'   ) > 0
#'
#'   if (has_ov) addHyperStableFraction(.) else {
#'     message("No overlap with hg38UltraStableProbes in this toy dataset; skipping.")
#'     .
#'   }
#' } %>%
#' qsea::getSampleTable() %>%
#' dplyr::select(sample_name, dplyr::any_of("hyperStableEdgar")) %>%
#' head()
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
    dplyr::filter(CpG_density >= !!minDensity) %>%
    dplyr::select(tidyselect::matches("beta"))

  overp8Fraction <- hyperStableBetaTable %>%
    replace(., is.na(.), 0) %>%
    {. >= minBeta} %>%
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
#' Return a tidy summary of the most relevant QC fields from a `qseaSet`
#' `sampleTable`. If core library metrics are missing, they are first added via
#' [addLibraryInformation()].
#'
#' @param qseaSet `qseaSet`.  
#'   Input object from which to extract QC fields.  
#'   **Default:** none (must be supplied).
#'
#' @details
#' If the `sampleTable` lacks library metrics (e.g., `valid_fragments`), the
#' function calls [addLibraryInformation()] to populate them. The output keeps
#' `sample_name` and any columns whose names match the regex
#' `"valid_fragment|relH|hyperStable|ichorTumo"` (if present), then orders rows
#' by `sample_name`.
#'
#' @return
#' A tibble with one row per sample containing:
#' - `sample_name`
#' - any available QC columns matching `valid_fragment`, `relH`, `hyperStable`,
#'   or `ichorTumo`.  
#' Columns not present in the `sampleTable` are silently omitted.
#'
#' @seealso
#' [qsea::getSampleTable()], [addLibraryInformation()], [addHyperStableFraction()]
#'
#' @examples
#' \donttest{
#' data("exampleTumourNormal", package = "mesa")
#'
#' # QC summary (adds library info if missing), ordered by sample_name
#' exampleTumourNormal %>%
#'   getSampleQCSummary() %>%
#'   head()
#' }
#'
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
