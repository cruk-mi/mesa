utils::globalVariables(c("chr", "seqnames", "start","end",".","annotation", "value","library_factor",
                         "ROI_start","ROI_end","value","window_start","window_end","window","sample_name",
                         "width","group","relH","CpG_density","strand","state","gene_name","map","CNV",
                         "EMSEMBL","SYMBOL","GENENAME","geneChr","geneStart","geneEnd","geneLength","n",
                         "adjPval","group1","group2","sample1","sample2",".up","hyperStableFractionp8","hyperStableFractionp9","hyperStableEdgar",
                         "isProperPair","isUnmappedQuery","isSupplementaryAlignment","isDuplicate","hasUnmappedMate",
                       "isNotPassingQualityControls","rname","pos","isize","MQ","mapq","isFirstMateRead","isPaired",
                       "cigar", ".rowID", "feature","annoShort","type",".comparison",".ext",".value","total_fragments",
                       "isSecondaryAlignment","ROI_ID","ID","ENSEMBL","deltaBeta",
                       "map_hg38_1000kb", "map_hg38_10kb", "map_hg38_50kb", "map_hg38_500kb",
                       "gc_hg38_1000kb","gc_hg38_50kb","gc_hg38_500kb","gc_hg38_10kb","score",
                       "log2FC","group1new","sample2new","counts","tumour","rowIndex","nUp","nDown","landscape","shortAnno","nOverCutoff",
                       "afterOverBackNum","initialOverBackNum", "qname","inOut","nSign","nStrands",
                       "chromosome_name","start_position","end_position", "input_file",
                       "topVarNumInput", "windowSdName", "pcaName",
                       "fnValue", "group2new","gap"))

#' This function takes a genomic ranges object or data frame, and lifts over from hg19 to hg38.
#'
#' @param grOrDf Either a genomic ranges object or a data frame (with seqnames column), with or without chr in the seqnames
#' @return A genomic ranges object object with coordinates in hg38.
#' @export
#' @seealso \code{rtracklayer::liftOver}
#'
#' @examples
#' \donttest{
#' # Example 1: liftover from a data.frame (mix of 'chr' and no 'chr')
#' df <- data.frame(
#'   seqnames = c("1", "chr2"),
#'   start    = c(100000, 200000),
#'   end      = c(100100, 200100)
#' )
#' gr_hg38 <- liftOverHg19(df)
#' gr_hg38
#' GenomeInfoDb::genome(gr_hg38)  # should show "hg38"
#'
#' # Example 2: GRanges input using the same coordinates as above
#' #' # Note: some intervals may not map; that's expected.
#' gr <- GenomicRanges::GRanges(df$seqnames, IRanges::IRanges(df$start, df$end))
#' liftOverHg19(gr)
#' }
liftOverHg19 <- function(grOrDf){

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "Package \"rtracklayer\" must be installed to use this function.",
      call. = FALSE
    )
  }

  regions <- grOrDf %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = ifelse(stringr::str_detect(seqnames,"chr"), seqnames, paste0("chr", seqnames))) %>%
    plyranges::as_granges()

  message("Performing coordinate liftover from hg19 to hg38")
  liftover_result <- rtracklayer::liftOver(regions, mesa::hg19ToHg38.over.chain)
  GenomeInfoDb::genome(liftover_result) <- "hg38"

  return(unlist(liftover_result))
}

.getMesaGenome <- function() {
  getOption("mesa_genome")
}

.getMesaMart <- function() {
  getOption("mesa_mart")
}

.getMesaTxDb <- function() {
  getOption("mesa_TxDb")
}

.getMesaAnnoDb <- function() {
  getOption("mesa_annoDb")
}

#' This function sets a default genome for the annotateWindows function to use
#'
#' Currently supporting "hg38" and "GRCh38", this sets the TxDb and annoDb options to appropriate files.
#'
#' @param genome String corresponding to a genome.
#' @return Invisibly returns \code{TRUE} on successful completion.
#' @export
#'
#' @examples
#' # Set the default genome to hg38
#' setMesaGenome("hg38")
#' getOption("mesa_genome")
#'
#' # Use the GRCh38 alias
#' setMesaGenome("GRCh38")
#' getOption("mesa_genome")
#'
#' # Reset/remove the default genome
#' setMesaGenome(NULL)
#' is.null(getOption("mesa_genome"))
setMesaGenome <- function(genome){
  options("mesa_genome" = genome)
  return(invisible(TRUE))
}

#' This function sets a default TxDb for the annotateWindows function to use.
#'
#' For instance "TxDb.Hsapiens.UCSC.hg38.knownGene" for human, or "TxDb.Mmusculus.UCSC.mm10.knownGene" for mouse.
#'
#' @param TxDb A TxDb package (as a string, or the object itself)
#' @return Invisibly returns \code{TRUE} on successful completion.
#' @export
#' @seealso The ChIPseeker function  used by annotateWindows is annotatePeak, \link[ChIPseeker]{annotatePeak}
#'
#' @examples
#' # Set default TxDb for human
#' setMesaTxDb("TxDb.Hsapiens.UCSC.hg38.knownGene")
#'
#' # Remove the default TxDb
#' setMesaTxDb(NULL)
setMesaTxDb <- function(TxDb){
  
  if(is.null(TxDb)){
    options("mesa_TxDb" = NULL)
    return(invisible(TRUE))
  }
  
  if (!requireNamespace(TxDb, quietly = TRUE)) {
    stop(
      glue::glue("Package {TxDb} must exist and be installed to set as default TxDb. Please install or correct and run again."),
      call. = FALSE
    )
  }

  options("mesa_TxDb" = TxDb)

  return(invisible(TRUE))
}

#' This function sets a default annoDb for the annotateWindows function to use.
#'
#' For instance "org.Hs.eg.db" for human, or "org.Mm.eg.db" for mouse.
#'
#' @param annoDb An annoDb package (as a string, or the object itself)
#' @return Invisibly returns \code{TRUE} on successful completion.
#' @export
#' @seealso The ChIPseeker function  used by annotateWindows is annotatePeak, \link[ChIPseeker]{annotatePeak}
#' 
#' @examples
#' # Set default annoDb for human (only if the package is installed)
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   setMesaAnnoDb("org.Hs.eg.db")
#'   getOption("mesa_annoDb")
#' }
#'
#' # Set default annoDb for mouse (only if the package is installed)
#' if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
#'   setMesaAnnoDb("org.Mm.eg.db")
#'   getOption("mesa_annoDb")
#' }
#'
#' # Unset the default annoDb
#' setMesaAnnoDb(NULL)
#' is.null(getOption("mesa_annoDb"))
setMesaAnnoDb <- function(annoDb){
  
  if(is.null(annoDb)){
    options("mesa_annoDb" = NULL)
    return(invisible(TRUE))
  }
  
  if (!requireNamespace(annoDb, quietly = TRUE)) {
    stop(
      glue::glue("Package {annoDb} must exist and be installed to set as default annoDb. Please install or correct and run again."),
      call. = FALSE
    )
  }
  
  options("mesa_annoDb" = annoDb)
  
  return(invisible(TRUE))
}

#' Manage mesa parallelisation (set & query)
#'
#' `setMesaParallel()` enables/disables parallelisation (and may register a backend);
#' `getMesaParallel()` returns the current setting.
#'
#' @param nCores How many cores to set. Assumes you want to use MulticoreParam, for another architecture use the useParallel argument and specify manually.
#' @param useParallel Boolean denoting whether or not to use parallelisation for various functions in mesa.
#' @param verbose Boolean to determine whether to print messages or not.
#' 
#' @return
#' - `setMesaParallel()` — logical scalar (TRUE/FALSE), returned **invisibly**; primarily used for its side effect
#'   of setting the `"mesa_parallel"` option and (optionally) registering a backend.
#' - `getMesaParallel()` — logical scalar indicating whether parallelisation is currently enabled; when
#'   `verbose = TRUE`, also prints a status message.
#'
#' @examples
#' # Turn parallelisation OFF (serial evaluation)
#' setMesaParallel(useParallel = FALSE, verbose = FALSE)
#' getMesaParallel()  # FALSE
#'
#' # --- Unix/mac (MulticoreParam): enable and use 2 cores ---
#' \donttest{
#' if (.Platform$OS.type != "windows") {
#'   old <- BiocParallel::bpparam()                         # save current backend
#'   setMesaParallel(nCores = 2, verbose = FALSE)           # registers MulticoreParam(2)
#'   BiocParallel::bpworkers()                               # e.g. 2
#'   getMesaParallel()                                       # TRUE
#'   BiocParallel::register(old)                             # restore previous backend
#'   setMesaParallel(useParallel = FALSE, verbose = FALSE)   # back to serial
#' }
#' }
#'
#' # --- Windows: register SnowParam yourself, then enable ---
#' \donttest{
#' if (.Platform$OS.type == "windows") {
#'   old <- BiocParallel::bpparam()
#'   BiocParallel::register(BiocParallel::SnowParam(workers = 2))
#'   setMesaParallel(useParallel = TRUE, verbose = FALSE)
#'   BiocParallel::bpworkers()                               # e.g. 2
#'   getMesaParallel()                                       # TRUE
#'   BiocParallel::register(old)
#'   setMesaParallel(useParallel = FALSE, verbose = FALSE)
#' }
#' }
#' @export
setMesaParallel <- function(nCores = NULL, useParallel = FALSE, verbose = TRUE){

  if(!is.null(nCores) && nCores > 1) {
    useParallel = TRUE
    BiocParallel::register(BiocParallel::MulticoreParam(workers = nCores))
  }

  options("mesa_parallel" = useParallel)

  if(useParallel & verbose){
    message(glue::glue("Parallelisation turned on for the functions in the mesa package, currently using {BiocParallel::bpworkers()} cores.
                       Control the number of cores by calling BiocParallel::register."))
  } else if (verbose) {
    message("Parallelisation turned off for all functions in the mesa package.")
  }

  return(invisible(useParallel))
}

#' @describeIn setMesaParallel Get whether mesa is using parallelisation.
#' @param verbose Boolean to determine whether to print messages or not.
#' @export
getMesaParallel <- function(verbose = FALSE) {

  if(verbose) {
    if (is.null(getOption("mesa_parallel"))) {
      message("Parallelisation not set, using serial evaluation. Call setMesaParallel to set for all mesa functions.")
    } else if (getOption("mesa_parallel")) {
      message(glue::glue("Using parallelisation, over {BiocParallel::bpworkers()} cores."))
    } else {
      message("Using serial evaluation")
    }
  }

  if (is.null(getOption("mesa_parallel"))) {
    return(FALSE)
  } else {
    return(getOption("mesa_parallel"))
  }

}

expect_no_error <- function(object) {
  testthat::expect_error({{ object }}, NA)
}
