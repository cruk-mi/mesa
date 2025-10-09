#' Drop fragment detail columns from regions
#'
#' Remove region metadata columns whose names contain `"avgFragment"` (e.g.,
#' `avgFragmentLength`, `avgFragmentMAPQ`) from the `regions` slot of a `qseaSet`,
#' if present.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object to modify.  
#'   **Default:** none (must be supplied).
#'
#' @details
#' This operation only affects the **metadata columns** (`mcols`) of the regions.
#' Genomic coordinates and the number/order of regions are preserved. If no
#' matching columns are found, the function is a no-op.
#'
#' @return
#' A `qseaSet` with its `regions` metadata updated (columns matching
#' `"avgFragment"` removed if present).
#'
#' @seealso
#' [qsea::getRegions()], [S4Vectors::mcols()], [dplyr::select()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda", package = "mesa") != "") {
#'   data(exampleTumourNormal, package = "mesa")
#'
#'   # Add toy fragment columns via plyranges, then drop them (pipe-only)
#'   exampleTumourNormal %>%
#'     { .@regions <- .@regions %>%
#'         plyranges::mutate(avgFragmentLength = 150L,
#'                           avgFragmentMAPQ   = 60L); . } %>%
#'     dropAvgFragDetails() %>%
#'     qsea::getRegions() %>%
#'     GenomicRanges::mcols() %>%
#'     colnames() %>%
#'     grep("^avgFragment", ., value = TRUE)   # character(0) if removed
#' }
#' }
#'
#' @export
dropAvgFragDetails <- function(qseaSet) {
  qseaSet@regions <- qseaSet@regions %>%
    dplyr::select(-tidyselect::matches("avgFragment"))

  return(qseaSet)
}


#' Convert a makeTable-like data frame to GRanges (UCSC style)
#'
#' Coerce a table of windows to a \linkS4class{GRanges}. Two input layouts are
#' supported:
#' - `chr`, `window_start`, `window_end` (as in `qsea::makeTable()` output), or
#' - `seqnames`, `start`, `end`.
#'
#' If chromosome labels lack the `"chr"` prefix, it is added automatically so
#' that `seqnames` are in UCSC style (e.g., `"chr1"`).
#'
#' @param dataTable `data.frame`.  
#'   Table of windows to convert. Must contain either the trio
#'   `chr`/`window_start`/`window_end` **or** `seqnames`/`start`/`end`.  
#'   **Default:** none (must be supplied).
#'
#' @details
#' - When given `chr`/`window_start`/`window_end`, the columns are renamed to
#'   `seqnames`/`start`/`end` before coercion.  
#' - When given `seqnames`/`start`/`end`, `seqnames` are first coerced to
#'   character and normalised to include the `"chr"` prefix.  
#' - The returned `GRanges` does **not** set a genome tag; set it yourself if
#'   needed (e.g., `GenomeInfoDb::genome(gr) <- "hg38"`).
#'
#' @return A `GRanges` with `seqnames` in UCSC style. The genome is unset.
#'
#' @seealso
#' [qsea::makeTable()], [plyranges::as_granges()],
#' [GenomeInfoDb::seqlevelsStyle()]
#'
#' @examples
#' # From makeTable-like columns
#' data.frame(chr = c("1","2"),
#'            window_start = c(100, 500),
#'            window_end   = c(200, 600)) %>%
#'   qseaTableToChrGRanges()
#'
#' # From seqnames/start/end; adds 'chr' if missing
#' data.frame(seqnames = c("chr3","4"),
#'            start = c(1000, 2000),
#'            end   = c(1100, 2100)) %>%
#'   qseaTableToChrGRanges()
#'
#' @export
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


#' Get window names from a qseaSet or ranges/table
#'
#' Return character labels of the form `"seqnames:start-end"` for each genomic
#' window/region.
#'
#' @param x `qseaSet` **or** `GRanges` **or** `data.frame`.  
#'   If a data frame, it must be coercible to `GRanges` (accepted columns include
#'   `seqnames/start/end`, or `chr/start/end`, or `chr/window_start/window_end`).  
#'   **Default:** none (must be supplied).
#'
#' @details
#' If `x` is a `qseaSet`, regions are taken from `qsea::getRegions(x)`. Otherwise
#' the input is coerced to `GRanges` (see [asValidGranges()]) and labels are
#' constructed as `"seqnames:start-end"`.
#'
#' @return
#' `character()` vector of window labels, one per region.
#'
#' @seealso
#' [qsea::getRegions()], [asValidGranges()]
#'
#' @examples
#' # From a GRanges (no intermediate objects)
#' GenomicRanges::GRanges(c("chr1","chr2"),
#'                        IRanges::IRanges(c(10,20), c(15,30))) %>%
#'   getWindowNames()
#'
#' # From a data.frame (no intermediate objects)
#' data.frame(seqnames = c("chr1","chr2"),
#'            start    = c(100,200),
#'            end      = c(120,220)) %>%
#'   getWindowNames()
#'
#' \donttest{
#' # From a qseaSet shipped with mesa (if available)
#' data(exampleTumourNormal, package = "mesa")
#' exampleTumourNormal %>% getWindowNames() %>% head()
#' }
#'
#' @export
getWindowNames <- function(x) {

  if(is.qseaSet(x)){
    return(x %>%
      qsea::getRegions() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(window = paste0(seqnames,":",start,"-",end)) %>%
      dplyr::pull(window))
  }

  return(asValidGranges(x) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(window = paste0(seqnames,":",start,"-",end)) %>%
    dplyr::pull(window))

  stop("Unknown data type!")
}


#' Infer pattern names from region density columns
#'
#' Scan the region metadata of a `qseaSet` for columns ending in `"_density"`
#' (e.g., `"CpG_density"`) and return the base pattern names (e.g., `"CpG"`).
#'
#' @param qseaSet `qseaSet`.  
#'   Object whose regions (via `qsea::getRegions()`) may contain `*_density`
#'   columns produced by `qsea::addPatternDensity()`.  
#'   **Default:** none (must be supplied).
#'
#' @details
#' This function inspects `mcols(qsea::getRegions(qseaSet))`, selects column
#' names matching the regex `"_density$"`, and strips that suffix. It does not
#' compute densities or modify the object. If no matching columns exist, a
#' zero-length character vector is returned.
#'
#' @return
#' `character()` vector of pattern names discovered (length 0 if none).
#'
#' @seealso
#' [qsea::addPatternDensity()], [qsea::getRegions()], [S4Vectors::mcols()]
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#' # Returns character(0) if no *_density columns are present
#' exampleTumourNormal %>% getPattern()
#' }
#'
#' @export
getPattern <- function(qseaSet) {
  qseaSet %>%
    qsea::getRegions() %>%
    GenomicRanges::mcols() %>%
    colnames() %>%
    stringr::str_subset("_density$") %>%
    stringr::str_remove("_density$") %>%
    return()
}


#' Extract the regions (windows) used in a qseaSet
#'
#' Convenience wrapper for [qsea::getRegions()], returning the genomic windows.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object.  
#'   **Default:** none (must be supplied).
#'
#' @return
#' A [GenomicRanges::GRanges] of windows.
#'
#' @seealso
#' [qsea::getRegions()], [getWindowNames()]
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#' exampleTumourNormal %>% getWindows() %>% head()
#' }
#'
#' @export
getWindows <- function(qseaSet) {
  qseaSet %>%
    qsea::getRegions() %>%
    return()
}


#' Remove almost-empty columns from a data frame
#'
#' Drop columns whose proportion of `NA` values exceeds a threshold.
#' Inspired by `janitor::remove_empty_cols()` but with a custom `prop` cutoff.
#'
#' @param dat `data.frame`.  
#'   Table to filter.  
#'   **Default:** none (must be supplied).
#'
#' @param prop `numeric(1)`.  
#'   Keep columns with `NA` proportion `<= prop` (range `[0, 1]`).  
#'   **Default:** none (must be supplied).
#'
#' @details
#' Column order and the number of rows are preserved. A message is emitted via
#' an internal janitor helper indicating which columns were removed.
#'
#' @return
#' A `data.frame` with columns failing the threshold removed.
#'
#' @seealso
#' [janitor::remove_empty()]
#'
#' @examples
#' \donttest{
#' data.frame(a = c(1, NA, 3),
#'            b = c(NA, NA, NA),
#'            c = c(1, 2, NA)) %>%
#'   remove_almost_empty_cols(prop = 0.5)  # drops column b
#' }
remove_almost_empty_cols <- function(dat, prop)  {
  mask_keep <- colSums(is.na(dat)) <=  prop*(nrow(dat))
  janitor:::remove_message(dat = dat, mask_keep = mask_keep, which = "cols", reason = "almost empty")
  return(dat[,mask_keep, drop = FALSE])
}


#' Skip slow checks when configured
#'
#' Test helper that skips long-running checks when the option
#' `options(skip_long_checks = TRUE)` is set.
#'
#' @return
#' Invisibly returns `TRUE` if **not** skipping; otherwise calls
#' `testthat::skip()` to skip the test.
#'
#' @examples
#' \donttest{
#' old <- options(skip_long_checks = TRUE)
#' skip_long_checks()  # will skip if option is TRUE
#' options(old)
#' }
#' @keywords internal
#' @noRd
skip_long_checks <- function() {
  if (!identical(options("skip_long_checks"), TRUE)) {
    return(invisible(TRUE))
  }

  testthat::skip("Slow checks skipped when options(skip_long_checks = TRUE) has been set")
}


#' Coerce common tabular inputs to GRanges
#'
#' Accept a `GRanges` or a data frame containing window coordinates and
#' return a valid [GenomicRanges::GRanges]. Supported data-frame schemas:
#' - `seqnames`, `start`, `end`
#' - `chr`, `start`, `end` (renamed to `seqnames`)
#' - `chr`, `window_start`, `window_end` (renamed to `seqnames/start/end`)
#'
#' @param object `GRanges` **or** `data.frame` with window coordinates.  
#'   **Default:** none (must be supplied).
#'
#' @return
#' A [GenomicRanges::GRanges] built from `object`.  
#' Errors if no supported schema is found.
#'
#' @seealso
#' [plyranges::as_granges()], [qseaTableToChrGRanges()]
#'
#' @examples
#' # From GRanges (passes through)
#' GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200)) %>%
#'   asValidGranges()
#'
#' # From data.frame (seqnames/start/end)
#' data.frame(seqnames = "chr2", start = 10, end = 20) %>%
#'   asValidGranges()
#'
#' # From data.frame (chr/window_start/window_end)
#' data.frame(chr = "1", window_start = 1000, window_end = 1100) %>%
#'   asValidGranges()
#' @export
asValidGranges <- function(object){

  if("GRanges" %in% class(object)){
    return(object)
  }

  if(is.data.frame(object)){
    if(length(intersect(colnames(object),c("seqnames","start","end"))) == 3){
      return(object %>% plyranges::as_granges())
    } else  if(length(intersect(colnames(object),c("chr","start","end"))) == 3){
      return(object %>% dplyr::rename(seqnames = chr) %>% plyranges::as_granges() )
    } else  if(length(intersect(colnames(object),c("chr","window_start","window_end"))) == 3){
      return(object %>%
               dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
               plyranges::as_granges() )
    } else {
      stop("Data frame can not be coerced to a GRanges object, requires seqnames, start, end columns.")
    }
  }

  stop("Object can not be coerced to a GRanges object")

}

#' This function captures a printed form of the object x, for use in error messages
#' Taken from https://stackoverflow.com/a/26083626 by Richie Cotton
#' @param x An object to capture
print_and_capture <- function(x) {paste(utils::capture.output(print(x)), collapse = "\n") }