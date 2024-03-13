#' Summarise the fraction of hyperstably methylated regions with clear methylation (hg38 only)

#' This function calculates and summarises the proportion of hyperstable regions that are present in each sample, 
#'   and adds this as a column to the sampleTable. Uses a set of windows shown to be methylated across a large number of 
#'   human methylation arrays over a wide range of tissues and diseases by [Edgar et al (2014)](https://pubmed.ncbi.nlm.nih.gov/25493099/)
#' @param qseaSet A qseaSet object
#' @export

addHyperStableFraction <- function(qseaSet){

  if(qsea:::getGenome(exampleTumourNormal) != "BSgenome.Hsapiens.NCBI.GRCh38") {
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

#' This function returns the most important columns of a qseaSet with respect to quality controls
#' @param qseaSet A qseaSet object
#' @export
#' @examples
#' exampleTumourNormal %>% getSampleQCSummary()
#' 
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
