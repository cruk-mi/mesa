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
                       "fnValue", "group2new","gap","size"))

#' Lift over genomic ranges from hg19 to hg38
#'
#' Convert genomic intervals from UCSC **hg19** to **hg38/GRCh38** using
#' a preloaded chain object.
#'
#' @param grOrDf `GRanges` or `data.frame`.  
#'   Input intervals. Data frames must contain at least `seqnames`, `start`,
#'   and `end`. `seqnames` may include or omit the `"chr"` prefix; it is added
#'   automatically if missing.  
#'   **Default:** none (must be supplied).
#'
#' @return
#' A \linkS4class{GRanges} in hg38 coordinates (with
#' `GenomeInfoDb::genome(x) <- "hg38"`). Intervals that cannot be mapped by the
#' chain are dropped (i.e., may return fewer ranges than provided).
#'
#' @details
#' Internally, the input is coerced to a `GRanges` (adding `"chr"` if needed),
#' then lifted over via \pkg{rtracklayer} with the chain stored in
#' `mesa::hg19ToHg38.over.chain`. Unmappable or split mappings are handled by
#' `rtracklayer::liftOver()`, and unmapped intervals are discarded when the
#' result is unlisted. The returned object is tagged as `"hg38"`.
#'
#' @seealso
#' [rtracklayer::liftOver()], [GenomicRanges::GRanges], [GenomeInfoDb::genome]
#'
#' @examples
#' # Data.frame input; mix of 'chr' and no 'chr' is fine
#' data.frame(
#'   seqnames = c("1", "chr2"),
#'   start    = c(100000, 200000),
#'   end      = c(100100, 200100)
#' ) %>% 
#'   liftOverHg19() %>% 
#'   { GenomeInfoDb::genome(.) }
#'
#' # GRanges input (some intervals may not map and can be dropped)
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100000, 100300))
#' liftOverHg19(gr)
#'
#' @export
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


#' Set default genome for downstream annotation helpers
#'
#' Record a preferred genome string (e.g., `"hg38"` / `"GRCh38"`) in the session
#' options so that functions like [annotateWindows()] can use it when `genome`
#' is not explicitly supplied.
#'
#' @param genome `character(1)` or `NULL`.  
#'   Genome identifier to store (typically `"hg38"` or `"GRCh38"`).  
#'   Use `NULL` to clear the setting.  
#'   **Default:** none (must be supplied).
#'
#' @details
#' This function sets the session option `options(mesa_genome = <genome>)`.
#' It **does not** modify TxDb or OrgDb options directly; set those via
#' [setMesaTxDb()] and [setMesaAnnoDb()] if you want global defaults for
#' transcript and gene annotation packages. Helpers such as [annotateWindows()]
#' consult `getOption("mesa_genome")` when `genome` is missing.
#'
#' @return
#' Invisibly returns `TRUE` on success.
#'
#' @seealso
#' [annotateWindows()], [setMesaTxDb()], [setMesaAnnoDb()]
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
#'
#' @export
setMesaGenome <- function(genome){
  options("mesa_genome" = genome)
  return(invisible(TRUE))
}


#' Set default TxDb for downstream annotation helpers
#'
#' Record a preferred TxDb to use when functions such as [annotateWindows()]
#' are called without an explicit `TxDb`.
#'
#' @param TxDb `character(1)`, TxDb **object**, or `NULL`.  
#'   Either the name of an installed TxDb package (e.g.,
#'   `"TxDb.Hsapiens.UCSC.hg38.knownGene"`), the TxDb object itself, or `NULL`
#'   to clear the setting. If a string is supplied, the package must be installed.  
#'   **Default:** none (must be supplied).
#'
#' @details
#' This sets the session option `options(mesa_TxDb = TxDb)`. When a string is
#' stored, helpers like [annotateWindows()] will load the TxDb object on demand
#' (e.g., `"TxDb.Hsapiens.UCSC.hg38.knownGene"`). Use [setMesaGenome()] and
#' [setMesaAnnoDb()] to configure complementary defaults for genome and OrgDb.
#'
#' @return
#' Invisibly returns `TRUE` on success.
#'
#' @seealso
#' [annotateWindows()], [setMesaGenome()], [setMesaAnnoDb()],
#' [ChIPseeker::annotatePeak()]
#'
#' @examples
#' # Set default TxDb for human (GRCh38/hg38)
#' setMesaTxDb("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' getOption("mesa_TxDb")
#'
#' # Clear the default TxDb
#' setMesaTxDb(NULL)
#' getOption("mesa_TxDb")
#'
#' # (Alternatively) you can store the object itself if already loaded:
#' # setMesaTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' @export
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


#' Set default OrgDb (annoDb) for downstream annotation helpers
#'
#' Record a preferred **OrgDb** package to use when functions such as
#' [annotateWindows()] are called without an explicit `annoDb`.
#'
#' @param annoDb `character(1)` or `NULL`.  
#'   Name of an installed OrgDb package (e.g., `"org.Hs.eg.db"` or
#'   `"org.Mm.eg.db"`), or `NULL` to clear the setting. The package must be
#'   installed if a string is supplied.  
#'   **Default:** none (must be supplied).
#'
#' @details
#' This sets the session option `options(mesa_annoDb = <annoDb>)`. When present,
#' helpers like [annotateWindows()] will use this OrgDb by default (particularly
#' for GRCh38 workflows) unless an explicit `annoDb` is provided. Use
#' [setMesaGenome()] and [setMesaTxDb()] to configure complementary defaults for
#' the genome and TxDb, respectively.
#'
#' @return
#' Invisibly returns `TRUE` on success.
#'
#' @seealso
#' [annotateWindows()], [setMesaGenome()], [setMesaTxDb()],
#' [ChIPseeker::annotatePeak()]
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
#'
#' @export
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
#' `setMesaParallel()` enables/disables parallelisation (and may register a
#' backend); `getMesaParallel()` returns the current setting.
#'
#' @param nCores `integer(1)` or `NULL`.  
#'   If `> 1` on Unix/macOS, registers `BiocParallel::MulticoreParam(workers = nCores)`
#'   and enables parallel mode. Ignored on Windows (see **Details**).  
#'   **Default:** `NULL`.
#'
#' @param useParallel `logical(1)`.  
#'   Explicitly turn mesa parallelisation on/off (sets the `"mesa_parallel"`
#'   option). If `nCores > 1`, this is treated as `TRUE`.  
#'   **Default:** `FALSE`.
#'
#' @param verbose `logical(1)`.  
#'   Print status messages.  
#'   **Default:** `TRUE`.
#'
#' @details
#' - This function stores the flag in `options(mesa_parallel = <TRUE/FALSE>)`.
#' - On Unix/macOS, supplying `nCores > 1` will **register** a
#'   `MulticoreParam` backend via `BiocParallel::register()` and enable
#'   parallelisation.
#' - On **Windows**, forked backends are unavailable; you should manually
#'   register a compatible backend (e.g., `SnowParam`) and then call
#'   `setMesaParallel(useParallel = TRUE)`.
#' - You can always inspect the active backend with `BiocParallel::bpparam()`
#'   and the worker count with `BiocParallel::bpworkers()`.
#'
#' @return
#' - `setMesaParallel()` — logical scalar (`TRUE`/`FALSE`), returned **invisibly**;
#'   primarily used for its side effects (setting the option and possibly
#'   registering a backend).
#' - `getMesaParallel()` — logical scalar indicating whether mesa parallelisation
#'   is currently enabled; with `verbose = TRUE`, also prints a status message.
#'
#' @seealso
#' [BiocParallel::register()], [BiocParallel::MulticoreParam()],
#' [BiocParallel::SnowParam()], [BiocParallel::bpworkers()]
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
#'
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
