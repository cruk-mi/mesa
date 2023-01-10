#' This function removes the avgFragmentLength and avgFragmentMAPQ columns from the regions of the qseaSet if present.
#' @param qseaSet A qseaSet object to drop the columns from the regions.
#'
dropAvgFragDetails <- function(qseaSet) {
  qseaSet@regions <- qseaSet@regions %>%
    dplyr::select(-tidyselect::matches("avgFragment"))

  return(qseaSet)
}

#' Convert a table from makeTable into a GRanges object
#'
#' This function adds the column names to be able to make the makeTable output into a hg38 GRanges object.
#' Basically it converts chr to seqnames, window_start to start, window_end to end, then turns it into a GRanges.
#'
#' @param dataTable The data frame, generally from qsea::makeTable. Can also be a GRanges object.
#' @return A GenomicRanges object

qseaTableToChrGRanges <- function(dataTable) {

  if ("window_start" %in% colnames(dataTable) ) {
    outGRanges <- dataTable %>%
      tibble::as_tibble() %>%
      dplyr::mutate(chr = as.character(chr)) %>%
      dplyr::mutate(seqnames = ifelse(stringr::str_detect(chr,"chr"),chr,paste0("chr", chr)), start = window_start, end = window_end) %>%
      plyranges::as_granges()

  } else if ("start" %in% colnames(dataTable)) {
    outGRanges <- dataTable %>%
      tibble::as_tibble() %>%
      dplyr::mutate(seqnames = as.character(seqnames)) %>%
      dplyr::mutate(seqnames = ifelse(stringr::str_detect(seqnames,"chr"),seqnames,paste0("chr", seqnames))) %>%
      plyranges::as_granges()
  } else {stop("Missing start or window_start column.")}

  return(outGRanges)
}

#' This function returns a vector of window names from a qseaSet or a data frame.
#' @param x A data frame or qseaSet to return the window names of
#' @export
#'
getWindowNames <- function(x) {

  if(is.qseaSet(x)){
    x %>%
      qsea::getRegions() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(window = paste0(seqnames,":",start,"-",end)) %>%
      dplyr::pull(window) %>%
      return()
  }

  if(is.data.frame(x) ){
    #TODO check that seqnames, start, end are present.
    x %>%
      tibble::as_tibble() %>%
      dplyr::mutate(window = paste0(seqnames,":",start,"-",end)) %>%
      dplyr::pull(window) %>%
      return()
  }
  #TODO method for GRanges object.
  else {stop("Unknown data type!")}
}

#' This function returns the name of the pattern added by addPatternDensity, e.g. "CpG"
#' @param qseaSet A qseaSet with pattern details added by addPatternDensity
#' @export
#' @examples
#' getPattern(exampleTumourNormal)
#'
getPattern <- function(qseaSet) {
  qseaSet %>%
    qsea::getRegions() %>%
    GenomicRanges::mcols() %>%
    colnames() %>%
    stringr::str_subset("_density$") %>%
    stringr::str_remove("_density$") %>%
    return()
}

#' This function removes columns that are mostly NA (with some proportion)
#' Based off the janitor::remove_empty_cols function.
#' @param dat A data frame to filter almost empty columns from
#' @param prop A proportion of NAs in each column to allow
#'
remove_almost_empty_cols <- function(dat, prop)  {
  mask_keep <- colSums(is.na(dat)) <=  prop*(nrow(dat))
  janitor:::remove_message(dat = dat, mask_keep = mask_keep, which = "cols", reason = "almost empty")
  return(dat[,mask_keep, drop = FALSE])
}

#' This function skips long running tests if options(skip_long_checks = TRUE) has been set
skip_long_checks <- function() {
  if (!identical(options("skip_long_checks"), TRUE)) {
    return(invisible(TRUE))
  }

  testthat::skip("Slow checks skipped when options(skip_long_checks = TRUE) has been set")
}
