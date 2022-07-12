#' This function calculates and summarises the proportion of hyperstable regions that are present in each sample, and adds to columns of the sampleTable
#' @param qseaSet A qseaSet object
#' @export
#'
addHyperStableFraction <- function(qseaSet){

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


#' This function returns the most important columns of a qseaSet with respect to quality controls
#' @param qseaSet A qseaSet object
#' @export
#'
getSampleQCSummary <- function(qseaSet){

  if (!("valid_fragments" %in% colnames(qsea::getSampleTable(qseaSet)))) {
    qseaSet <- qseaSet %>%
      addLibraryInformation()
  }

  if (!("hyperStableFractionp8" %in% colnames(qsea::getSampleTable(qseaSet)))) {
    qseaSet <- qseaSet %>%
      addHyperStableFraction()
  }

  qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::select(sample_name, tidyselect::matches("valid_fragment|relH|hyperStable|ichorTumo")) %>%
    dplyr::arrange(sample_name) %>%
    tibble::as_tibble() %>%
    return()

}

