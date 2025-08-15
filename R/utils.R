#' Drop fragment detail columns from regions
#'
#' Removes region metadata columns whose names match `"avgFragment"` (e.g.,
#' `avgFragmentLength`, `avgFragmentMAPQ`) from the `regions` slot of a `qseaSet`,
#' if present.
#'
#' @param qseaSet A `qseaSet` object to modify.
#'
#' @return The input `qseaSet` with its `regions` metadata updated (columns matching
#' `"avgFragment"` removed if present).
#'
#' @seealso [qsea::getRegions()], [plyranges::as_granges()], [dplyr::select()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package = "mesa")
#'   # Add toy fragment columns, then drop them
#'   qs <- exampleTumourNormal
#'   r  <- qsea::getRegions(qs)
#'   GenomicRanges::mcols(r)$avgFragmentLength <- 150L
#'   GenomicRanges::mcols(r)$avgFragmentMAPQ   <- 60L
#'   qs@regions <- r
#'   qs2 <- dropAvgFragDetails(qs)
#'   setdiff(colnames(GenomicRanges::mcols(qs@regions)),
#'           colnames(GenomicRanges::mcols(qs2@regions)))
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
#' Coerces a table of windows to `GRanges`. Two input formats are supported:
#' - `chr`, `window_start`, `window_end` (like `qsea::makeTable()` output), or
#' - `seqnames`, `start`, `end`.
#' If chromosome labels lack `"chr"`, the prefix is added.
#'
#' @param dataTable A data frame (typically from `qsea::makeTable()`), or a `GRanges`.
#'
#' @return A `GRanges` with `seqnames` in UCSC style (`"chr1"`, …). The genome tag
#' is not set; set it yourself if needed (e.g., `GenomeInfoDb::genome(gr) <- "hg38"`).
#'
#' @seealso [qsea::makeTable()], [plyranges::as_granges()], [GenomeInfoDb::seqlevelsStyle()]
#'
#' @examples
#' # From makeTable-like columns
#' df1 <- data.frame(chr = c("1","2"), window_start = c(100, 500), window_end = c(200, 600))
#' gr1 <- qseaTableToChrGRanges(df1)
#' gr1
#'
#' # From seqnames/start/end; adds chr prefix if missing
#' df2 <- data.frame(seqnames = c("chr3","4"), start = c(1000, 2000), end = c(1100, 2100))
#' gr2 <- qseaTableToChrGRanges(df2)
#' gr2
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
#' Returns character labels of the form `"seqnames:start-end"` for each window.
#'
#' @param x A `qseaSet`, a `GRanges`, or a data frame coercible to `GRanges`
#' (must contain `seqnames/chr`, `start`, `end`).
#'
#' @return A character vector of window labels, one per region.
#'
#' @seealso [qsea::getRegions()], [asValidGranges()]
#'
#' @examples
#' # From a GRanges
#' gr <- GenomicRanges::GRanges(c("chr1","chr2"), IRanges::IRanges(c(10,20), c(15,30)))
#' getWindowNames(gr)
#'
#' # From a data frame
#' df <- data.frame(seqnames = c("chr1","chr2"), start = c(100,200), end = c(120,220))
#' getWindowNames(df)
#'
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   head(getWindowNames(exampleTumourNormal))
#' }
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

#' Infer the pattern name used by addPatternDensity
#'
#' Scans the region metadata of a `qseaSet` for columns ending with `"_density"`
#' and returns the base pattern name (e.g., `"CpG"` for `"CpG_density"`).
#'
#' @param qseaSet A `qseaSet` with pattern density columns (e.g., added by `addPatternDensity()`).
#'
#' @return A character vector of pattern names found (possibly length 0 if none).
#'
#' @seealso [qsea::getRegions()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   # Returns character(0) if no *_density columns present
#'   getPattern(exampleTumourNormal)
#' }
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
#' Convenience wrapper for `qsea::getRegions()`, returning the genomic windows.
#'
#' @param qseaSet A `qseaSet`.
#'
#' @return A `GRanges` of windows.
#'
#' @seealso [qsea::getRegions()], [getWindowNames()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   getWindows(exampleTumourNormal)
#' }
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
#' Drops columns whose proportion of `NA` values exceeds a threshold.
#' Based on `janitor::remove_empty_cols()`, but with a custom `prop` cutoff.
#'
#' @param dat A data frame.
#' @param prop Numeric in `[0, 1]`. Keep columns with `NA` proportion `<= prop`.
#'
#' @return A data frame with columns failing the threshold removed. The original
#' order and row count are preserved.
#'
#' @seealso [janitor::remove_empty()]
#'
#' @examples
#' df <- data.frame(a = c(1, NA, 3),
#'                  b = c(NA, NA, NA),
#'                  c = c(1, 2, NA))
#' remove_almost_empty_cols(df, prop = 0.5)  # drops column b
#' 
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
#' @return Invisibly returns `TRUE` if not skipping; otherwise calls `testthat::skip()`.
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
#' Checks if `object` is a `GRanges` or a data frame containing window coordinates,
#' and returns a valid `GRanges`. Supported data frame schemas:
#' - `seqnames`, `start`, `end`
#' - `chr`, `start`, `end` (renamed to `seqnames`)
#' - `chr`, `window_start`, `window_end` (renamed to `seqnames/start/end`)
#'
#' @param object A `GRanges` or a data frame with window coordinates.
#'
#' @return A `GRanges` with coordinates taken from `object`. Errors if no supported
#' schema is found.
#'
#' @seealso [plyranges::as_granges()], [qseaTableToChrGRanges()]
#'
#' @examples
#' # From GRanges (passes through)
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
#' asValidGranges(gr)
#'
#' # From data frame (seqnames/start/end)
#' df1 <- data.frame(seqnames = "chr2", start = 10, end = 20)
#' asValidGranges(df1)
#'
#' # From data frame (chr/window_start/window_end)
#' df2 <- data.frame(chr = "1", window_start = 1000, window_end = 1100)
#' asValidGranges(df2)
#' 
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
