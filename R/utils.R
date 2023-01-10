#' This function removes the avgFragmentLength and avgFragmentMAPQ columns from the regions of the qseaSet
#' @param qseaSet A qseaSet object to drop the columns from the regions.
#' @export
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

#' This function returns the window names from a data frame with DMR windows
#' @param dataTable A data frame to return the window names from
#' @export
#'
getWindowNames <- function(dataTable) {
  dataTable %>%
    dplyr::mutate(window = paste0(seqnames,":",start,"-",end)) %>%
    dplyr::pull(window) %>%
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
