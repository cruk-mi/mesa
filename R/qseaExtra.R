# Internal: memoise getCGPositions so that it is cached as it takes a long time
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  getCGPositions <<- memoise::memoise(getCGPositions)
}

#' Check whether an object is a qseaSet
#'
#' Predicate helper used in input validation and branching logic.
#'
#' @param x `ANY`.  
#'   Object to test.
#'
#' @return `logical(1)`: `TRUE` if `x` inherits from class `"qseaSet"`, otherwise `FALSE`.
#'
#' @seealso
#' [qsea::createQseaSet()], [qsea::getSampleTable()], [base::inherits()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # TRUE for qseaSet objects
#' exampleTumourNormal %>% is.qseaSet()
#'
#' # FALSE for non-qseaSet objects
#' iris %>% is.qseaSet()
#'
#' @export
is.qseaSet <- function(x){
  return(inherits(x,"qseaSet"))
  }

#' Set or get an Ensembl/BioMart handle on a qseaSet
#'
#' Store and retrieve a “mart” handle inside a `qseaSet` for downstream
#' annotation tasks (e.g., resolving gene coordinates).
#'
#' @name setMart
#'
#' @param object `qseaSet`.  
#'   The object to modify or query.
#'
#' @param mart `biomaRt::Mart`, `character(1)`, or `ANY`.  
#'   Connection object returned by **biomaRt** (recommended), or a string label
#'   you interpret elsewhere (e.g., `"ENSEMBL_110"`). Stored verbatim in
#'   `object@parameters$mart`.  
#'   **Default:** none (must be supplied to `setMart()`).
#'   
#' @param ... Additional arguments passed to methods or to
#'   [biomaRt::useMart()]. Not used in the default `"qseaSet"` method.
#'
#' @return
#' - `setMart()` returns the updated `qseaSet` (with `@parameters$mart` set).  
#' - `getMart()` returns the stored value (whatever was set), or `NULL` if absent.
#'
#' @details
#' These are S4 generics with methods for class `"qseaSet"`. The value is not
#' validated or dereferenced here; functions that consume it (e.g., gene
#' annotation helpers) should perform any needed checks.
#'
#' @section Methods:
#' \describe{
#'   \item{`setMart(object, mart)`}{Set the handle on a `qseaSet`.}
#'   \item{`getMart(object)`}{Retrieve the stored handle from a `qseaSet`.}
#' }
#'
#' @seealso
#' [annotateWindows()], [setMesaGenome()], [setMesaTxDb()], [setMesaAnnoDb()],
#' [biomaRt::useMart()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Store a label (or use a real biomaRt::Mart) and read it back
#' exampleTumourNormal %>%
#'    setMart("ENSEMBL_110") %>%
#'    getMart()
NULL

#' @rdname setMart
#' @export
setGeneric('setMart', function(object, ...) standardGeneric('setMart'))

#' @rdname setMart
#' @export
setMethod('setMart', 'qseaSet', function(object, mart){
  object@parameters$mart <- mart
  object
})

#' @rdname setMart
#' @export
setGeneric('getMart', function(object, ...) standardGeneric('getMart'))

#' @rdname setMart
#' @export
setMethod('getMart', 'qseaSet', function(object) object@parameters$mart)

#' Annotate genomic windows using ChIPseeker (with optional CpG/FANTOM context)
#'
#' Uses [ChIPseeker::annotatePeak()] to assign transcript-relative annotations
#' (promoter, exon, intron, intergenic, etc.) to genomic windows in `dataTable`.
#' If a genome is specified or globally set via [setMesaGenome()], defaults and
#' convenience datasets are applied for GRCh38. CpG island context and FANTOM
#' enhancer overlaps can be added if the corresponding ranges are provided (or
#' left `NULL` to use mesa defaults for GRCh38).
#'
#' @param dataTable `data.frame`/`tibble` (coercible to `GRanges`) or `GRanges`.  
#'   Windows to annotate. Data frames must include `seqnames`, `start`, and `end`
#'   (or columns convertible by [qseaTableToChrGRanges()]).
#'
#' @param genome `character(1)` or `NULL`.  
#'   Guides annotation defaults. Currently supports `"hg38"`/`"GRCh38"`.  
#'   **Default:** value set by [setMesaGenome()] (via internal `.getMesaGenome()`), or `NULL`.
#'
#' @param TxDb TxDb object or `character(1)`.  
#'   Either an unquoted TxDb object or a string like
#'   `"TxDb.Hsapiens.UCSC.hg38.knownGene"`. If character, it is resolved at
#'   runtime.  
#'   **Default:** value set by [setMesaTxDb()] (via `.getMesaTxDb()`), or for
#'   GRCh38/hg38 when `NULL`, use
#'   `TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene`
#'   if installed.
#'
#' @param annoDb `character(1)` or `NULL`.  
#'   OrgDb package name (e.g., `"org.Hs.eg.db"`).  
#'   **Default:** value set by [setMesaAnnoDb()] (via `.getMesaAnnoDb()`), or for
#'   GRCh38/hg38 when `NULL`, use `"org.Hs.eg.db"` if installed.
#'
#' @param CpGislandsGR `GRanges` or `NULL`.  
#'   CpG island regions for island/shore/shelf context.  
#'   **Default:** `NULL` (for GRCh38/hg38, uses `mesa::hg38CpGIslands`).
#'
#' @param FantomRegionsGR `GRanges` or `NULL`.  
#'   FANTOM enhancer regions for overlap counts.  
#'   **Default:** `NULL` (for GRCh38/hg38, uses `mesa::FantomRegions`).
#'
#' @details
#' * **Conversion:** Non-`GRanges` inputs are converted via [qseaTableToChrGRanges()].  
#' * **Genome-aware defaults (GRCh38/hg38):** Missing `TxDb`/`annoDb`/contexts are
#'   filled using the packages noted above when installed.  
#' * **Seqlevels style:** Output `seqnames` have the `"chr"` prefix removed
#'   (e.g., `"chr1"` → `"1"`).  
#' * **ChIPseeker call:** Uses `tssRegion = c(-2000, 500)`, `level = "transcript"`,
#'   `overlap = "all"`, `verbose = FALSE`.  
#' * **Output columns:** Includes `annotation` (full label), `shortAnno` (without the
#'   trailing parenthetical), and—when contexts are provided—`landscape`
#'   (Island/Shore/Shelf/Open Sea) and `inFantom` (overlap count with FANTOM).
#'
#' @return A tibble with the input windows augmented by ChIPseeker annotations and,
#'   when available, CpG island landscape and Fantom overlap.
#'
#' @seealso
#' [ChIPseeker::annotatePeak()], [setMesaGenome()], [setMesaTxDb()], [setMesaAnnoDb()],
#' [qseaTableToChrGRanges()], [liftOverHg19()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Derive some regions (e.g., DMRs) then annotate using GRCh38 defaults
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "tumour", contrasts = "first") %>%
#'   annotateWindows(genome = "hg38")
#'
#' # Or specify TxDb/annoDb explicitly
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "tumour", contrasts = "first") %>%
#'   annotateWindows(
#'     TxDb   = "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'     annoDb = "org.Hs.eg.db"
#'   )
#'
#' # Mouse example (mm10): supply mouse TxDb and OrgDb
#' data(exampleMouse, package = "mesa")
#' exampleMouse %>%
#'   getRegions() %>%
#'   annotateWindows(
#'     TxDb   = "TxDb.Mmusculus.UCSC.mm10.knownGene",
#'     annoDb = "org.Mm.eg.db"
#'   )
#'
#' # You can also set defaults globally using setMesaTxDb and setMesaAnnoDb:
#'   setMesaGenome("hg38"); 
#'   setMesaTxDb("TxDb.Hsapiens.UCSC.hg38.knownGene"); 
#'   setMesaAnnoDb("org.Hs.eg.db")
#'
#' @export
annotateWindows <- function(dataTable, genome = .getMesaGenome(), TxDb = .getMesaTxDb(), 
                            annoDb = .getMesaAnnoDb(), CpGislandsGR = NULL,
                            FantomRegionsGR = NULL) {

  if(!is.null(TxDb) & is.character(TxDb)){
    TxDb <- eval(parse(text=paste0(TxDb,"::", TxDb)))
  }
  
  if(is.null(TxDb) & is.null(genome)) {
    stop("Please specify a TxDb or genome, this can be set globally using setMesaTxDb and/or setMesaGenome")
  }
  
  if(is.null(genome)){
    genome <- ""
  }
  
  if(genome %in% c("hg38","GRCh38") && is.null(TxDb)) {

    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      stop(
        "Package \"TxDb.Hsapiens.UCSC.hg38.knownGene\" must be installed to use this function. Please install and run again.",
        call. = FALSE
      )
    }
    TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }

  if(genome  %in% c("hg38","GRCh38") && is.null(annoDb)) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop(
        "Package \"org.Hs.eg.db\" must be installed to use this function. Please install and run again.",
        call. = FALSE
      )
    }
    annoDb <- "org.Hs.eg.db"
    }

  if(is.null(annoDb) && is.null(genome)) {
    stop("Please specify a annoDb or genome, this can be set globally using setMesaannoDb and/or setMesaGenome")
  }
  
  if(genome  %in% c("hg38","GRCh38") && is.null(CpGislandsGR)) { CpGislandsGR <- mesa::hg38CpGIslands }

  if(genome  %in% c("hg38","GRCh38") && is.null(FantomRegionsGR)) { FantomRegionsGR <- mesa::FantomRegions %>% plyranges::as_granges()}

  if(methods::is(dataTable,"GRanges")) {
    GRangesObject <- dataTable
  } else{
    GRangesObject <- dataTable %>%
      qseaTableToChrGRanges()
  }

  chipseekerData <- GRangesObject %>%
    ChIPseeker::annotatePeak(tssRegion = c(-2000, 500),
                             level = "transcript", # changed from gene to transcript to stop it outputting some genes as being >10Mb long
                             TxDb = TxDb,
                             annoDb = annoDb,
                             overlap = "all",
                             verbose = FALSE)

  grAnno <- chipseekerData@anno %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = stringr::str_remove(seqnames, "chr")) %>%
    plyranges::as_granges()


  if(!is.null(CpGislandsGR)){
  grAnno <- grAnno %>%
    plyranges::mutate(nIslands = plyranges::count_overlaps(., CpGislandsGR),
           nShore = plyranges::count_overlaps(., plyranges::flank_left(CpGislandsGR, width = 2000)) +
                    plyranges::count_overlaps(., plyranges::flank_right(CpGislandsGR, width = 2000)),
           nShelf = plyranges::count_overlaps(., plyranges::shift_left(plyranges::flank_left(CpGislandsGR, width = 2000), 2000)) +
                    plyranges::count_overlaps(., plyranges::shift_right(plyranges::flank_right(CpGislandsGR, width = 2000), 2000))) %>%
    dplyr::mutate(landscape = dplyr::case_when(nIslands > 0 ~ "Island",
                                               nShore > 0 ~ "Shore",
                                               nShelf > 0 ~ "Shelf",
                                               TRUE ~ "Open Sea"))
  }

  if(!is.null(FantomRegionsGR)){
   grAnno <- grAnno %>%
    plyranges::mutate(inFantom = plyranges::count_overlaps(., FantomRegionsGR))
  }

  dfAnno <- grAnno %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(shortAnno = stringr::str_remove(annotation," \\(.*")) %>%
    dplyr::relocate(dplyr::ends_with("means"), .after = dplyr::last_col()) %>%
    dplyr::select(-width, -strand)

  return(dfAnno)
}


#' Subset windows in a qseaSet by signal across samples
#'
#' Filter genomic windows based on a summary function applied across selected
#' samples. Typical uses include keeping windows with median NRPM above a
#' threshold, or retaining windows where the minimum beta is below 0.5.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing windows and signal.
#'
#' @param fn `function`.  
#'   Summary function applied **row-wise** across the selected samples (e.g.,
#'   `median`, `min`, `max`). If you need NA handling, provide it explicitly,
#'   e.g. `function(x) median(x, na.rm = TRUE)`.  
#'   **Default:** none (must be supplied).
#'
#' @param threshold `numeric(1)`.  
#'   Threshold applied to the summary value returned by `fn`.  
#'   **Default:** none (must be supplied).
#'
#' @param aboveThreshold `logical(1)`.  
#'   Keep windows with values `>= threshold` when `TRUE`; keep windows with values
#'   `< threshold` when `FALSE`.  
#'   **Default:** none (must be supplied).
#'
#' @param samples `character()` or `NULL`.  
#'   Sample names to include, or a single string used as a pattern to match
#'   sample names. If `NULL`, use all samples.  
#'   **Default:** `NULL`.
#'
#' @param normMethod `character(1)`.  
#'   Normalisation/measure to evaluate. One of `"nrpm"`, `"beta"`, or `"counts"`.  
#'   **Default:** `"nrpm"`.
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, aggregate replicates by the `group` column in the sample table
#'   and apply `fn` to group means; if `FALSE`, use individual samples.  
#'   **Default:** `FALSE`.
#'
#' @details
#' The function builds a window × sample matrix for the chosen `normMethod`,
#' optionally restricts to samples via `samples` (exact names or pattern match),
#' computes `fn` per row, and filters by `threshold` according to
#' `aboveThreshold`. Supply a summary function that accepts a numeric vector and
#' returns a scalar; include `na.rm = TRUE` inside `fn` if needed.
#'
#' @return A filtered `qseaSet` containing only the selected windows.
#'
#' @seealso
#' [getDataTable()], [filterByOverlaps()]
#'
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Keep windows with median NRPM > 1 across all samples
#' exampleTumourNormal %>%
#'   subsetWindowsBySignal(
#'     fn = median,
#'     threshold = 1,
#'     aboveThreshold = TRUE
#'   )
#'
#' # Keep windows where the MIN beta < 0.5 in any sample
#' exampleTumourNormal %>%
#'   subsetWindowsBySignal(
#'     fn = min,
#'     threshold = 0.5,
#'     aboveThreshold = FALSE,
#'     normMethod = "beta"
#'   )
#'
#' # Restrict to samples whose names contain "Lung"
#' exampleTumourNormal %>%
#'   subsetWindowsBySignal(
#'     fn = median,
#'     threshold = 1,
#'     aboveThreshold = TRUE,
#'     samples = "Lung"
#'   )
#'
#' # Use group means instead of individual samples
#' exampleTumourNormal %>%
#'   subsetWindowsBySignal(
#'     fn = median,
#'     threshold = 1,
#'     aboveThreshold = TRUE,
#'     useGroupMeans = TRUE
#'   )
#'
#' @export
subsetWindowsBySignal <- function(qseaSet, fn, threshold, aboveThreshold, samples = NULL, normMethod = "nrpm", useGroupMeans = FALSE){

  fnName <- as.character(substitute(fn, env = environment()))

  if (!useGroupMeans) {
    qseaSamples <- qsea::getSampleNames(qseaSet)
    groupString <- ""
  } else {
    qseaSamples <- names(getSampleGroups2(qseaSet))
    groupString <- " group"
  }

  samplesNotInQset <- setdiff(samples, qseaSamples)

  if (length(samples) == 1 & length(samplesNotInQset) > 0 & is.character(samples)) {
    sampleNameString <- samples
    samples <- stringr::str_subset(qseaSamples, samples)
    message(glue::glue("Considering {length(samples)} sample{groupString}s containing \"{sampleNameString}\" in the name."))
  } else if (length(samplesNotInQset) > 0 ) {
    stop(glue::glue("Sample{groupString}(s) {paste0(samplesNotInQset, collapse = ', ')} not present in the qseaSet!"))
  }

  if (is.null(samples)) {
    samples <- qseaSamples
    message(glue::glue("Considering all {length(samples)} sample{groupString}s."))
  }

  if (length(samples) == 0) {
    stop("No samples selected.")
  }

  if (!length(normMethod(normMethod)) == 1) {
    stop(glue::glue("normMethod should be a single valid option for qsea::normMethod"))
  }

  if (!useGroupMeans) {
    dataTable <- getDataTable(qseaSet %>% dplyr::filter(sample_name %in% !!samples), normMethod = normMethod, useGroupMeans = useGroupMeans)
  } else {
    dataTable <- getDataTable(qseaSet %>% dplyr::filter(group %in% !!samples), normMethod = normMethod, useGroupMeans = useGroupMeans)
  }

  dataTable <- withCallingHandlers(
    dataTable %>%
      dplyr::mutate(
        fnValue = apply(dplyr::pick(tidyselect::all_of(samples)), 1, fn, na.rm = TRUE)
      ),
    warning = function(w) {
      if (grepl("no non-missing arguments to .*; returning -Inf", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )

  if (aboveThreshold) {
    dataTable <- dataTable %>% dplyr::filter(fnValue >= !!threshold)
    keepString <- "above (or equal to)"
  } else {
    dataTable <- dataTable %>% dplyr::filter(fnValue < !!threshold)
    keepString <- "below"
  }

  message(glue::glue("Keeping {nrow(dataTable)} windows with {fnName} {keepString} {threshold} over {length(samples)} sample{groupString}s."))

  qseaSet <- filterByOverlaps(qseaSet, dataTable)
  return(qseaSet)

}


#' Subset windows above Poisson background
#'
#' Identify genomic windows with significantly more reads than expected under a
#' Poisson background model, given total reads per sample and the number of
#' genome windows. If `numWindows` is provided, then that value will be used as
#' the total original number of windows under consideration for the false 
#' discovery rate, otherwise it will use the number of windows currently in the
#' qseaSet.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing per-window read counts.
#'
#' @param keepAbove `logical(1)`.  
#'   If `TRUE`, **keep** windows above background; if `FALSE`, **remove** them.  
#'   **Default:** `FALSE`.
#'
#' @param samples `character()` or `NULL`.  
#'   Vector of sample names to test, or a single string used as a pattern matched
#'   against sample names. If `NULL`, all samples are tested.  
#'   **Default:** `NULL`.
#'
#' @param numWindows `integer(1)` or `NULL`.  
#'   Total number of windows in the genome used to compute the expected background.
#'   If `NULL`, it will use all the windows currently in the qseaSet.
#'   **Default:** `NULL`.
#'
#' @param FDRthres `numeric(1)`.  
#'   FDR threshold used to call windows significantly above background.  
#'   **Default:** `0.01`.
#'
#' @param numAbove `integer(1)`.  
#'   Minimum number of tested samples that must exceed background in a window
#'   to keep/drop it (depending on `keepAbove`).  
#'   **Default:** `1`.
#'
#' @details
#' For each tested sample, an expected count per window is derived from its total
#' reads and the supplied/estimated `numWindows`, and a Poisson test is used to
#' flag windows above background. This uses the `valid_fragments` column in the 
#' library information attached to the qseaSet, which is the number of fragments
#' used in the generation of the original qseaSet. Windows are retained or 
#' removed based on `keepAbove` and the requirement that at least `numAbove` 
#' samples are significant at `FDRthres`.
#'
#' @return A filtered `qseaSet` containing only windows passing the criteria.
#'
#' @seealso
#' [subsetWindowsBySignal()]
#'
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Keep windows above Poisson background in at least 2 samples whose names 
#' # contain "Lung", assuming there were 9 million original windows
#' exampleTumourNormal %>%
#'   subsetWindowsOverBackground(keepAbove = TRUE, samples = "Lung", 
#'                               numAbove = 2, numWindows = 9e6 )
#'
#' # Drop windows above background (retain only background-like windows) across 
#' # all samples, assuming there were 9 million original windows
#' exampleTumourNormal %>%
#'   subsetWindowsOverBackground(keepAbove = FALSE, numWindows = 9e6)
#' @importFrom rlang :=
#' @export
subsetWindowsOverBackground <- function(qseaSet, keepAbove = FALSE, 
                                        samples = NULL, numWindows = NULL,
                                        FDRthres = 0.01, numAbove = 1){
  
  samplesNotInSet <- setdiff(samples, qsea::getSampleNames(qseaSet))
  
  if (length(samples) == 1 & length(samplesNotInSet) > 0 & is.character(samples)) {
    sampleNameString <- samples
    samples <- stringr::str_subset(qsea::getSampleNames(qseaSet), samples)
    message(glue::glue("Filtering out {length(samples)} samples containing {sampleNameString} in the name."))
  } else if (length(samplesNotInSet) > 0 ) {
    stop(glue::glue("Sample {samplesNotInSet} not present in the qseaSet!

                    "))
  }
  
  if (is.null(samples)) {
    samples <- qsea::getSampleNames(qseaSet)
  }
  
  message(glue::glue("Removing windows with reads above background levels in {length(samples)} samples."))
  
  countMat <- qseaSet %>%
    qsea::getCounts()
  
  if (is.null(numWindows)) {
    numWindows <- nrow(countMat)
  }
  
  fdrMat <- purrr::map_dfc(samples,
                           function(x){
                             totalNumReads <- qseaSet@libraries$file_name[x, "valid_fragments"]
                             lambda <- totalNumReads/numWindows
                             pvals <- stats::ppois(countMat[,x] - 1, lambda, lower.tail = FALSE)
                             fdrvals <- stats::p.adjust(pvals, method = "fdr") %>%
                               tibble::enframe(name = "window") %>%
                               dplyr::rename(!!x := value) %>%
                               dplyr::select(-window) 
                           }
  )
  
  if (keepAbove) {
    indexToKeep <- fdrMat %>%
      {. <= FDRthres} %>%
      rowSums() %>%
      {. >= numAbove} %>%
      which()
  } else {
    indexToKeep <- fdrMat %>%
      {. <= FDRthres} %>%
      rowSums() %>%
      {. < numAbove} %>%
      which()
  }
  
  windowsToKeep <- qsea::getRegions(qseaSet)[indexToKeep]
  
  print(windowsToKeep)
  
  message(glue::glue("Removing {nrow(fdrMat) - length(windowsToKeep)} windows based on {length(samples)} samples, {length(windowsToKeep)} remaining"))
  
  return(filterByOverlaps(qseaSet, windowsToKeep))
  
}


#' Downsample reads in a qseaSet
#'
#' Randomly subsample reads to a fixed depth across all samples in a `qseaSet`.
#' Useful for equalising library sizes prior to comparison.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object whose per-window counts will be downsampled.  
#'   **Default:** none (must be supplied).
#'
#' @param nReads `integer(1)`.  
#'   Target number of reads to retain **per sample**. If a sample has fewer than
#'   `nReads`, it is left unchanged.  
#'   **Default:** none (must be supplied).
#'
#' @return
#' A `qseaSet` with downsampled counts. Library metadata
#' (`valid_fragments`, `offset`, `library_factor`) are updated accordingly.
#'
#' @details
#' Downsampling is performed independently per sample. To obtain reproducible
#' results across runs, call `set.seed()` **before** invoking `downSample()`.
#' Samples with total reads `< nReads` are not upsampled.
#'
#' @seealso
#' [getDataTable()], [subsetWindowsBySignal()], [mixSamples()]
#'
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' set.seed(1)
#' # Downsample to 100 fragments per sample, then check per-sample totals
#' exampleTumourNormal %>%
#'   downSample(100) %>%
#'   getCountTable() %>%
#'   dplyr::select(tidyselect::matches("_[NT]")) %>%
#'   colSums() %>%
#'   range()  # should be ~ c(100, 100) if all samples had >= 100 reads
#'
#' # Also inspect the updated library sizes recorded in the sample table
#' exampleTumourNormal %>%
#'   downSample(100) %>%
#'   qsea::getSampleTable() %>%
#'   dplyr::select(patient, type, gender)
#'
#' @export
downSample <- function(qseaSet, nReads){
  counts <- qseaSet@count_matrix

  if (min(colSums(counts)) < nReads) {
    stop(glue::glue("Number of reads requested is less than the minimum {min(colSums(counts))}."))
  }

  message(glue::glue("Downsampling all samples to {nReads} each"))

  newCounts <- purrr::map_dfc(colnames(counts), function(colname){
        vec <- counts[,colname]

        sample(rep(seq_along(vec), vec), replace = FALSE, size = nReads) %>%
          table() %>%
          tibble::enframe(name = "window") %>%
          dplyr::mutate(window = as.integer(window), value = as.integer(value)) %>%
          dplyr::left_join(tibble::tibble(window = seq_along(vec)),., copy = TRUE, by = "window") %>%
          dplyr::mutate(value = tidyr::replace_na(value,0)) %>%
          dplyr::select(value) %>%
          dplyr::rename(!!colname := value)

       }
     )

  qseaSet@count_matrix <- as.matrix(newCounts)
  qseaSet@libraries$file_name[,"valid_fragments"] <- rep(nReads, ncol(counts))
  qseaSet@libraries$file_name[,"offset"] <- rep(NA, ncol(counts))
  qseaSet@libraries$file_name[,"library_factor"] <- rep(NA, ncol(counts))

  return(qseaSet)
}


#' Convert qsea beta values to array-like probe matrix
#'
#' Build a probe × sample beta table by mapping `qseaSet` windows to array
#' probe coordinates. Currently supports Illumina **Infinium 450k** probes on
#' GRCh38; a custom `GRanges` of probes can also be supplied.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object providing per-window beta values.  
#'   **Default:** none (must be supplied).
#'
#' @param arrayDetails `character(1)` or `GRanges`.  
#'   Either a recognised keyword (currently `"Infinium450k"`) or a `GRanges`
#'   of probe loci with an `ID` metadata column of probe identifiers. Coordinates
#'   must match the `qseaSet` genome (e.g., GRCh38).  
#'   **Default:** `"Infinium450k"`.
#'
#' @return
#' A `data.frame` with one row per probe (including an `ID` column) and one
#' column per sample containing the corresponding beta values. Probes with no
#' matching/overlapping window receive `NA`.
#'
#' @details
#' Probe coordinates are matched to `qseaSet` windows (e.g., via overlap). When
#' multiple windows are relevant for a probe, the package's implementation is
#' used to derive a single beta value per probe. Ensure the probe genome build
#' matches the `qseaSet` (e.g., GRCh38) to avoid mismatches.
#'
#' @seealso
#' [getDataTable()], dataset `mesa::hg38_450kArrayGR`
#'
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Use the built-in keyword for 450k probes (GRCh38)
#' exampleTumourNormal %>%
#'   convertToArrayBetaTable(arrayDetails = "Infinium450k") %>%
#'   head()
#'
#' # Supply an explicit GRanges of 450k probes (must have metadata column 'ID')
#' exampleTumourNormal %>%
#'   convertToArrayBetaTable(arrayDetails = mesa::hg38_450kArrayGR) %>%
#'   head()
#'
#' @export
convertToArrayBetaTable <- function(qseaSet, arrayDetails = "Infinium450k") {

  if(is.character(arrayDetails)){
    if(arrayDetails == "Infinium450k" & qsea:::getGenome(qseaSet) == "BSgenome.Hsapiens.NCBI.GRCh38"){
      arrayObject <- mesa::hg38_450kArrayGR
    } else if(arrayDetails == "Infinium450k" & qsea:::getGenome(qseaSet) == "BSgenome.Hsapiens.UCSC.hg38"){
      arrayObject <- mesa::hg38_450kArrayGR %>% tibble::as_tibble() %>% dplyr::mutate(seqnames = paste0("chr",seqnames)) %>% plyranges::as_granges()
    }
    else {stop("Only Infinium450k implemented currently as a string.")}

  } else {
    arrayObject <- asValidGranges(arrayDetails)
  }

  if(!("ID" %in% colnames(GenomicRanges::mcols(arrayObject)))){
    stop("No ID column found in object")
  }

  qseaSet %>%
    qsea::makeTable(norm_methods = "beta", samples = qsea::getSampleNames(.), ROIs = arrayObject) %>%
    dplyr::select(-tidyselect::matches("ROI_start|ROI_end|ROI_chr")) %>%
    dplyr::rename(ID = ROI_ID) %>%
    dplyr::select(ID, tidyselect::matches("beta")) %>%
    dplyr::rename_with( ~ stringr::str_remove(., "_beta"), tidyselect::matches("_beta"))

}


#' Fraction of thresholded windows overlapping a set of regions
#'
#' For each sample, compute the proportion of **windows** with counts
#' `>= numCountsNeeded` that overlap `regionsToOverlap`, relative to all
#' windows with counts `>= numCountsNeeded` in that sample.
#'
#' @details
#' Operates on the window × sample count matrix from [qsea::getCounts()].
#' Counts are converted to a logical indicator (count `>= numCountsNeeded`)
#' and summed per sample:
#' - `initialOverBackNum`: windows meeting the threshold genome-wide.
#' - `afterOverBackNum`: windows meeting the threshold **within**
#'   `regionsToOverlap` (via [filterByOverlaps()]).
#' The reported `fraction` is `afterOverBackNum / initialOverBackNum`.
#'
#' **Note:** despite the function name, this computes a **window-based** fraction
#' using a count threshold, not a direct fraction of raw reads.
#'
#' @param qseaSet `qseaSet`.  
#'   Input qseaSet. **Default:** none (must be supplied).
#'
#' @param regionsToOverlap `GRanges` or `data.frame`.  
#'   Regions to consider for overlap. Data frames must be coercible to `GRanges`
#'   (e.g., have `seqnames`/`start`/`end` or `chr`/`window_start`/`window_end`).  
#'   **Default:** none (must be supplied).
#'
#' @param numCountsNeeded `integer(1)`.  
#'   Minimum reads per window for it to count toward the fraction.  
#'   **Default:** none (must be supplied).
#'
#' @return
#' A tibble with one row per sample containing:
#' - `sample_name`  
#' - `initialOverBackNum` — windows `>=` threshold genome-wide  
#' - `afterOverBackNum` — windows `>=` threshold within `regionsToOverlap`  
#' - `fraction` — `afterOverBackNum / initialOverBackNum`  
#' followed by columns from the sample table (left-joined by `sample_name`).
#'
#' @seealso
#' [qsea::getCounts()], [qsea::getSampleTable()], [filterByOverlaps()],
#' [hg38CpGIslands], [subsetWindowsBySignal()]
#'
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Using the shipped ultra-stable probes (GRCh38) as the region set
#' exampleTumourNormal %>%
#'   calculateFractionReadsInGRanges(
#'     regionsToOverlap = mesa::hg38CpGIslands,
#'     numCountsNeeded   = 5
#'   ) %>%
#'   dplyr::select(sample_name, initialOverBackNum, afterOverBackNum, fraction) 
#'
#' # Define a small GRanges subset from the object's own windows
#' gr <- exampleTumourNormal %>%
#'   qsea::getRegions() %>%
#'   head(n = 100)
#'   
#' calculateFractionReadsInGRanges(
#'       qseaSet           = exampleTumourNormal,
#'       regionsToOverlap  = gr,
#'       numCountsNeeded   = 3) %>%
#'   dplyr::select(sample_name, initialOverBackNum, afterOverBackNum, fraction) 
#'
#' @export
calculateFractionReadsInGRanges <- function(qseaSet, regionsToOverlap, numCountsNeeded) {
  initialReadTotals <- qseaSet %>%
    qsea::getCounts() %>%
    {. >= numCountsNeeded } %>%
    colSums()

  afterSubsetReadTotals <- qseaSet %>%
    filterByOverlaps(regionsToOverlap) %>%
    qsea::getCounts() %>%
    {. >= numCountsNeeded } %>%
    colSums()

  out <- tibble::tibble(sample_name = names(initialReadTotals),
         initialOverBackNum = initialReadTotals,
         afterOverBackNum = afterSubsetReadTotals,
         fraction = afterOverBackNum/initialOverBackNum) %>%
    dplyr::left_join(qsea::getSampleTable(addLibraryInformation(qseaSet)))
  
  return(out)

}


#' Remove normalisation suffix from column names
#'
#' Clean wide tables (e.g., from [getDataTable()]) by stripping the trailing
#' `_{normMethod}` or `_{normMethod}_means` suffixes from sample /
#' group-mean columns.
#'
#' @param dataTable `data.frame` or tibble.  
#'   Window × sample (or group) table whose column names may include a
#'   normalisation suffix (e.g., `"sampleA_beta"`, `"Tumour_beta_means"`).  
#'   **Default:** none (must be supplied).
#'
#' @param normMethod `character(1)`.  
#'   Normalisation/measure tag to remove from column names (e.g., `"beta"`,
#'   `"nrpm"`).  
#'   **Default:** none (must be supplied).
#'
#' @return
#' A `data.frame` like `dataTable` but with cleaned column names (suffixes
#' `_{normMethod}` and `_{normMethod}_means` removed where present).
#'
#' @details
#' Intended for post-processing tables where column names encode the measure,
#' e.g., `"sample_beta"`, `"group_beta_means"`. The function drops those suffixes
#' to yield bare sample/group names, facilitating downstream joins and plotting.
#'
#' @seealso
#' [getDataTable()], [selectQset()], [pull.qseaSet()]
#'
#' @family table-helpers
#'
#' @examples
#' # Basic cleaning for 'beta'
#' data.frame(A_beta = 1:3, B_beta_means = 4:6) %>%
#'   removeNormMethodSuffix("beta") %>%
#'   names()
#'
#' # Works similarly for other measures (e.g., 'nrpm')
#' data.frame(S1_nrpm = 1:2, Tumour_nrpm_means = 3:4) %>%
#'   removeNormMethodSuffix("nrpm") %>%
#'   names()
#'
#' @export
removeNormMethodSuffix <- function(dataTable, normMethod) {
  dplyr::rename_with(dataTable, ~ stringr::str_remove(.x, glue::glue("_{normMethod}(_means)?$")))
}


#' Count windows above a cutoff
#'
#' For a given set of genomic windows, count (per sample) how many windows
#' exceed a chosen threshold using the selected normalisation/measure.
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object containing methylation-enriched sequencing data.  
#'   **Default:** none (must be supplied).
#'
#' @param GRanges `GenomicRanges::GRanges`.  
#'   Windows to count over.  
#'   **Default:** none (must be supplied).
#'
#' @param samples `character()` or `NULL`.  
#'   Sample names to include, or a single string used as a pattern to match
#'   sample names. If `NULL`, all samples are used.  
#'   **Default:** `NULL`.
#'
#' @param cutoff `numeric(1)`.  
#'   Threshold; windows with value `>= cutoff` are counted.  
#'   **Default:** `0`.
#'
#' @param normMethod `character(1)`.  
#'   Measure to use. One of `"nrpm"` or `"beta"`.  
#'   **Default:** `"nrpm"`.
#'
#' @details
#' The function builds a window × sample matrix over `GRanges` for the chosen
#' `normMethod`, optionally restricts to `samples`, and counts, per sample, the
#' number of windows with value `>= cutoff`. It also reports how many windows
#' were evaluated (`totalWindowsUsed`). Results are joined with the sample table.
#'
#' @return
#' A tibble with one row per sample containing:
#' - `sample_name`  
#' - `numOverCutoff` — number of windows `>= cutoff`  
#' - `totalWindowsUsed` — number of windows evaluated  
#' followed by columns from the sample table (left-joined).
#'
#' @seealso
#' [getDataTable()], [summariseAcrossWindows()]
#'
#' @family window-summaries
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' #' # Count windows with NRPM >= 1 over the first 100 regions
#' countWindowsAboveCutoff(exampleTumourNormal, 
#'                                qsea::getRegions(exampleTumourNormal)[1:100], 
#'                                cutoff = 1, 
#'                                normMethod = "nrpm")
#'                                
#' @export
countWindowsAboveCutoff <- function(qseaSet, GRanges, samples = NULL,
                                   cutoff = 0, normMethod = "nrpm"){

  if (is.null(samples)) {
    samples <- qsea::getSampleNames(qseaSet)
  }

  reducedData <- qseaSet %>%
    filterByOverlaps(GRanges) %>%
    qsea::makeTable(norm_methods = normMethod, samples = samples) %>%
    dplyr::select(dplyr::matches(normMethod)) %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "_beta|_nrpm|_means")) %>%
    {. >= cutoff}

  reducedData %>%
    colSums() %>%
    tibble::enframe(name = "sample_name", value = "numOverCutoff") %>%
    dplyr::mutate(totalWindowsUsed = nrow(reducedData)) %>%
    dplyr::left_join(qsea::getSampleTable(qseaSet))

}


#' Make a wide sample-by-window table
#'
#' Create a wide table for ML/visualisation where **rows are samples** and
#' **columns are genomic windows**, using one normalisation/measure.
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object containing methylation-enriched sequencing data. 
#'   **Default:** none (must be supplied).
#'
#' @param normMethod `character(1)`.  
#'   Measure to extract. Common choices include `"nrpm"` or `"beta"`.  
#'   **Default:** `"nrpm"`.
#'
#' @param ... tidyselect specification or `character()`.  
#'   Optional columns from the sample table to append as additional features
#'   (e.g., `group, tumour` or `c("group","tumour")`).  
#'   **Default:** none.
#'
#' @details
#' Internally obtains a window × sample matrix for `normMethod`, transposes it
#' to sample × window, and (if requested) appends selected sample-table columns.
#' Window columns are typically named `"chr:start-end"`. If measure suffixes
#' are present in column names, they may be cleaned to bare sample/window names
#' (see [removeNormMethodSuffix()]).
#'
#' @return
#' A tibble with one row per sample and one column per genomic window, plus any
#' appended sample-level metadata.
#'
#' @seealso
#' [getDataTable()], [removeNormMethodSuffix()], [countWindowsAboveCutoff()]
#'
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # NRPM features (sample × first few windows) with a group column appended
#' exampleTumourNormal %>%
#'   makeTransposedTable(normMethod = "nrpm", group) %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Beta features; inspect dimensions
#' exampleTumourNormal %>%
#'   makeTransposedTable(normMethod = "beta") %>%
#'   dim()
#'
#' @export
makeTransposedTable <- function(qseaSet, normMethod = "nrpm", ...){

  #TODO: Rewrite to use getDataTable rather than makeTable.
  
  qseaSet %>%
    qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = normMethod) %>%
    dplyr::rename_with(~ stringr::str_replace_all(.x, "_nrpm$|_beta$|_counts$", ""))  %>%
    dplyr::select(-CpG_density) %>%
    tidyr::pivot_longer(-c(chr, window_start,window_end), names_to = "sample_name", values_to = "value") %>%
    dplyr::mutate(chr = ifelse(stringr::str_detect(chr,"chr"), chr, paste0("chr",chr))) %>%
    tidyr::unite(col = "window", chr, window_start, window_end) %>%
    tidyr::pivot_wider(names_from = window, values_from = value) %>%
    dplyr::left_join(qsea::getSampleTable(qseaSet) %>% dplyr::select(sample_name, ...)) %>%
    dplyr::relocate(tidyselect::matches("^chr"), .after = tidyselect::last_col())

}


#' Get counts per window
#'
#' Convenience wrapper around [getDataTable()] to extract **raw counts**.
#' Returns a window × sample (or group) table for downstream summaries/plots.
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'   **Default:** none (must be supplied).
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, average replicates by the `group` column in the sample table and
#'   return group means; if `FALSE`, return per-sample counts.  
#'   **Default:** `FALSE`.
#'
#' @param addMethodSuffix `logical(1)`.  
#'   If `TRUE`, keep the method suffix in column names (e.g., `Sample1_counts`
#'   or `Group_counts_means`); if `FALSE`, return bare sample/group names.  
#'   **Default:** `FALSE`.
#'
#' @param verbose `logical(1)`.  
#'   Print progress/messages.  
#'   **Default:** `TRUE`.
#'   
#' @return
#' A tibble/data.frame with one row per genomic window and columns for each
#' sample (or group), plus any window metadata included by [getDataTable()].
#'
#' @seealso
#' [getDataTable()], [getNRPMTable()], [getBetaTable()]
#'
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Per-sample counts (first rows/columns)
#' exampleTumourNormal %>%
#'   getCountTable() %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Group means instead of individual samples
#' exampleTumourNormal %>%
#'   getCountTable(useGroupMeans = TRUE) %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Keep the method suffix in column names
#' exampleTumourNormal %>%
#'   getCountTable(addMethodSuffix = TRUE) %>%
#'   names() %>%
#'   head()
#'
#' @rdname getCountTable
#' @export
getCountTable <- function(qseaSet, useGroupMeans = FALSE, addMethodSuffix = FALSE, verbose = TRUE){
  tab <- qseaSet %>% 
    getDataTable(normMethod = "counts", 
                 useGroupMeans = useGroupMeans, 
                 addMethodSuffix = addMethodSuffix,
                 verbose = verbose)
  
  return(tab)
}


#' Get NRPM per window
#'
#' Convenience wrapper around [getDataTable()] to extract **NRPM** values
#' (normalised reads per million). Returns a window × sample (or group) table
#' for downstream summaries/plots.
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'   **Default:** none (must be supplied).
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, average replicates by the `group` column in the sample table and
#'   return group means; if `FALSE`, return per-sample NRPM.  
#'   **Default:** `FALSE`.
#'
#' @param addMethodSuffix `logical(1)`.  
#'   If `TRUE`, keep the method suffix in column names (e.g., `Sample1_nrpm`
#'   or `Group_nrpm_means`); if `FALSE`, return bare sample/group names.  
#'   **Default:** `FALSE`.
#'
#' @param verbose `logical(1)`.  
#'   Print progress/messages.  
#'   **Default:** `TRUE`.
#'
#' @return
#' A tibble/data.frame with one row per genomic window and columns for each
#' sample (or group), plus any window metadata included by [getDataTable()].
#'
#' @seealso
#' [getDataTable()], [getCountTable()], [getBetaTable()]
#'
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Per-sample NRPM (first rows/columns)
#' exampleTumourNormal %>%
#'   getNRPMTable() %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Group means instead of individual samples
#' exampleTumourNormal %>%
#'   getNRPMTable(useGroupMeans = TRUE) %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Keep the method suffix in column names
#' exampleTumourNormal %>%
#'   getNRPMTable(addMethodSuffix = TRUE) %>%
#'   names() %>%
#'   head()
#'
#' @rdname getNRPMTable
#' @export
getNRPMTable <- function(qseaSet, useGroupMeans = FALSE, addMethodSuffix = FALSE, verbose = TRUE){
  tab <- qseaSet %>% 
    getDataTable(normMethod = "nrpm", 
                 useGroupMeans = useGroupMeans, 
                 addMethodSuffix = addMethodSuffix, 
                 verbose = verbose)
  
  return(tab)
}


#' Get beta per window
#'
#' Convenience wrapper around [getDataTable()] to extract **beta** values
#' (fraction methylated). Returns a window × sample (or group) table for
#' downstream summaries/plots.
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object containing methylation-enriched sequencing data. 
#'   **Default:** none (must be supplied).
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, average replicates by the `group` column in the sample table and
#'   return group means; if `FALSE`, return per-sample betas.  
#'   **Default:** `FALSE`.
#'
#' @param minEnrichment `integer(1)`.  
#'   Minimum reads per window required for a non-`NA` beta (below this threshold,
#'   qsea sets beta to `NA`).  
#'   **Default:** `3`.
#'
#' @param addMethodSuffix `logical(1)`.  
#'   If `TRUE`, keep the method suffix in column names (e.g., `Sample1_beta`
#'   or `Group_beta_means`); if `FALSE`, return bare sample/group names.  
#'   **Default:** `FALSE`.
#'
#' @param verbose `logical(1)`.  
#'   Print progress/messages.  
#'   **Default:** `TRUE`.
#'
#' @return
#' A tibble/data.frame with one row per genomic window and columns for each
#' sample (or group), plus any window metadata included by [getDataTable()].
#'
#' @seealso
#' [getDataTable()], [getCountTable()], [getNRPMTable()]
#'
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Per-sample beta (first rows/columns), with a minimum read threshold
#' exampleTumourNormal %>%
#'   getBetaTable(minEnrichment = 3) %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Group means instead of individual samples
#' exampleTumourNormal %>%
#'   getBetaTable(useGroupMeans = TRUE, minEnrichment = 3) %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Keep the method suffix in column names
#' exampleTumourNormal %>%
#'   getBetaTable(addMethodSuffix = TRUE) %>%
#'   names() %>%
#'   head()
#'
#' # Inspect how many NAs arise from the minEnrichment filter
#' exampleTumourNormal %>%
#'   getBetaTable(minEnrichment = 5) %>%
#'   dplyr::select(dplyr::where(is.numeric)) %>%
#'   is.na() %>%
#'   colSums() %>%
#'   head()
#'
#' @rdname getBetaTable
#' @export
getBetaTable <- function(qseaSet, useGroupMeans = FALSE, minEnrichment = 3, addMethodSuffix = FALSE, verbose = TRUE){
   tab <- qseaSet %>% 
     getDataTable(normMethod = "beta", 
                  useGroupMeans = useGroupMeans, 
                  minEnrichment = minEnrichment, 
                  addMethodSuffix = addMethodSuffix,
                  verbose = verbose)
  return(tab)
}


#' Summarise across windows per sample
#'
#' Apply a summary function (e.g., mean, median, sd) over genomic windows
#' per sample, optionally restricted to a set of regions and for one or more
#' normalisation methods.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object providing window-level signal.  
#'   **Default:** none (must be supplied).
#'
#' @param regionsToOverlap `GRanges`, `data.frame`, or `NULL`.  
#'   Regions used to restrict the windows before summarising. Data frames must be
#'   coercible to `GRanges` (e.g., have `seqnames`/`start`/`end` or
#'   `chr`/`window_start`/`window_end`). If `NULL`, use all windows in `qseaSet`.  
#'   **Default:** `NULL`.
#'
#' @param fn `function`.  
#'   Summary function applied **column-wise per sample** over the selected windows
#'   (e.g., `mean`, `median`, `sd`). Provide NA handling inside `fn` if needed,
#'   e.g., `function(x) mean(x, na.rm = TRUE)`.  
#'   **Default:** `mean`.
#'
#' @param suffix `character(1)`.  
#'   Optional string appended to output column names (e.g., a region label).  
#'   **Default:** `""`.
#'
#' @param addSampleTable `logical(1)`.  
#'   If `TRUE`, join sample-table columns to the result.  
#'   **Default:** `TRUE`.
#'
#' @param normMethod `character()`.  
#'   One or more measures to summarise. Typically a subset of `c("nrpm","beta")`.  
#'   **Default:** `c("nrpm","beta")`.
#'
#' @param naMethod `character(1)`.  
#'   How to treat missing values prior to/within summarisation. Common choices:  
#'   `"na.rm"` — call `fn` with `na.rm = TRUE` when supported;  
#'   `"drop"`  — drop windows with any `NA` across the selected methods before
#'   summarising.  
#'   **Default:** `"na.rm"`.
#'
#' @param minEnrichment `integer(1)`.  
#'   Minimum reads per window required for a non-`NA` beta (relevant when
#'   `normMethod` includes `"beta"`).  
#'   **Default:** `3`.
#'
#' @param fnName `character(1)` or `NULL`.  
#'   Name of `fn` used to construct output column names; if `NULL`, inferred
#'   from `fn`.  
#'   **Default:** `NULL`.
#'
#' @details
#' The function builds a window × sample matrix for each `normMethod`, optionally
#' restricts to `regionsToOverlap`, then applies `fn` over windows for each
#' sample to produce a single summary value per sample (per method). Missing-value
#' handling follows `naMethod`. When summarising betas, windows below
#' `minEnrichment` reads are `NA` prior to summarisation.
#'
#' @return
#' A tibble with one row per sample. For each requested `normMethod`, a summary
#' column is added (column naming typically reflects `<method>` and `fn` and may
#' include `suffix`). If `addSampleTable = TRUE`, sample-table columns are joined.
#'
#' @seealso
#' [addSummaryAcrossWindows()], [getDataTable()]
#'
#' @family window-summaries
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Mean NRPM and beta across all windows
#' exampleTumourNormal %>%
#'   summariseAcrossWindows(
#'     fn         = mean,
#'     normMethod = c("nrpm","beta")
#'   )
#'
#' # Maximum NRPM within a small genomic region
#' exampleTumourNormal %>%
#'   summariseAcrossWindows(
#'     regionsToOverlap = data.frame(seqnames = 7, start = 25002001, end = 25017900),
#'     fn         = max,
#'     normMethod = "nrpm",
#'     suffix     = "_chr7_slice"
#'   )
#'
#' # Median beta over all windows
#' exampleTumourNormal %>%
#'   summariseAcrossWindows(
#'     fn         = median,
#'     normMethod = "beta",
#'     minEnrichment = 3
#'   )
#'
#' @export
summariseAcrossWindows <- function(qseaSet,
                                   regionsToOverlap = NULL,
                                   fn = mean,
                                   addSampleTable = TRUE,
                                   normMethod = c("nrpm", "beta"),
                                   naMethod = "na.rm",
                                   minEnrichment = 3,
                                   suffix = "",
                                   fnName = NULL) {
    #TODO: Can we get multiple summary statistics in one go?

  if(is.null(fnName)) {
    fnName <- as.character(substitute(fn, env = environment()))
  }

    #if suffix doesn't start with "_" then add that to the string
    suffix <- ifelse(stringr::str_detect(suffix, "^_") | nchar(suffix) == 0 , suffix, paste0("_",suffix))

    if(is.null(regionsToOverlap)) {
      regionsToOverlap <- qsea::getRegions(qseaSet)
    }

    dataMat <- qseaSet %>%
      filterByOverlaps(regionsToOverlap) %>%
      getDataTable(normMethod = normMethod,
                   minEnrichment = minEnrichment,
                   addMethodSuffix = TRUE)

    if(naMethod == "drop"){

      message("Dropping rows with an NA value in any sample")

      dataMat <- dataMat %>%
        tidyr::drop_na()

    }

    if(naMethod == "na.rm"){
      message("Removing NA values on a per-sample basis")
    }

    map_out <- purrr::map(normMethod,
                   function(normType){
                     temp <- dataMat %>%
                       dplyr::select(dplyr::matches(paste0("_", normType,"$"))) %>%
                       dplyr::rename_with(~ stringr::str_remove_all(.x, paste0("_", normType, "$")))

                     out <- temp %>%
                       apply(2, fn, na.rm = TRUE) %>%
                       tibble::enframe(name = "sample_name", value = paste0(normType, "_", fnName, suffix)) %>%
                       dplyr::mutate(num_windows = temp %>% apply(2,function(x) !is.na(x)) %>% colSums()) %>%
                       dplyr::rename_with(~stringr::str_replace(.x, "num_windows", paste0(normType, "_num_windows", suffix)))

                   }
    ) %>%
      purrr::reduce(dplyr::full_join, by = "sample_name")

    if(addSampleTable) {
      map_out <- map_out %>%
        dplyr::left_join(qsea::getSampleTable(qseaSet), by = "sample_name")
    }

    return(map_out)

}


#' Append window summaries to the sample table
#'
#' Compute per-sample summaries across windows and left-join them to the
#' `qseaSet` sample table.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object providing window-level signal.  
#'   **Default:** none (must be supplied).
#'
#' @param regionsToOverlap `GRanges`, `data.frame`, or `NULL`.  
#'   Regions used to restrict the windows before summarising. Data frames must be
#'   coercible to `GRanges` (e.g., have `seqnames`/`start`/`end` or
#'   `chr`/`window_start`/`window_end`). If `NULL`, use all windows in `qseaSet`.  
#'   **Default:** `NULL`.
#'
#' @param fn `function`.  
#'   Summary function applied **per sample** over the selected windows (e.g.,
#'   `mean`, `median`, `sd`). Provide NA handling inside `fn` if needed, e.g.,
#'   `function(x) mean(x, na.rm = TRUE)`.  
#'   **Default:** `mean`.
#'
#' @param suffix `character(1)`.  
#'   Optional string appended to output column names (e.g., a region label) so
#'   you can call this repeatedly for different regions.  
#'   **Default:** `""`.
#'
#' @param normMethod `character()`.  
#'   One or more measures to summarise. Typically a subset of `c("nrpm","beta")`.  
#'   **Default:** `c("nrpm","beta")`.
#'
#' @param naMethod `character(1)`.  
#'   How to treat missing values prior to/within summarisation. Supported values:
#'   `"na.rm"` — call `fn` with `na.rm = TRUE` (when supported);  
#'   `"drop"`  — drop windows with any `NA` across selected methods;  
#'   `"impute"` — replace missing values via the package’s internal strategy
#'   before summarising.  
#'   **Default:** `"impute"`.
#'
#' @param minEnrichment `integer(1)`.  
#'   Minimum reads per window required for a non-`NA` beta (relevant when
#'   `normMethod` includes `"beta"`).  
#'   **Default:** `3`.
#'
#' @details
#' Internally, this calls [summariseAcrossWindows()] to compute per-sample
#' summaries for the requested `normMethod`(s), then left-joins those summaries
#' to the `sampleTable`. Output column names typically encode the method and the
#' summary function (and include `suffix` when provided).
#'
#' @return
#' A `qseaSet` with new summary columns appended to its `sampleTable`.
#'
#' @seealso
#' [summariseAcrossWindows()], [getDataTable()]
#'
#' @family window-summaries
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Add mean NRPM and beta across all windows 
#' exampleTumourNormal %>%
#'   addSummaryAcrossWindows(
#'     fn         = mean,
#'     normMethod = c("nrpm","beta")
#'   ) %>%
#'   getSampleTable()
#'
#' # Maximum NRPM across a small genomic slice
#' exampleTumourNormal %>%
#'   addSummaryAcrossWindows(
#'     regionsToOverlap = data.frame(seqnames = 7, start = 25002001, end = 25017900),
#'     fn         = max,
#'     normMethod = "nrpm"
#'   ) %>%
#'   getSampleTable()
#'
#' # Add summaries repeatedly with different regions using a suffix
#' exampleTumourNormal %>%
#'   addSummaryAcrossWindows(
#'     regionsToOverlap = data.frame(seqnames = 7, start = 25002001, end = 25017900),
#'     fn         = median,
#'     normMethod = "nrpm",
#'     suffix     = "subset"
#'   ) %>%
#'   addSummaryAcrossWindows(
#'     fn         = median,
#'     normMethod = "nrpm",
#'     suffix     = "all"
#'   ) %>%
#'   getSampleTable()
#'
#' @export
addSummaryAcrossWindows <- function(qseaSet,
                                    regionsToOverlap = NULL,
                                    fn = mean,
                                    suffix = "",
                                    normMethod = c("nrpm", "beta"),
                                    naMethod = "impute",
                                    minEnrichment = 3) {

  #need to catch function name when called like this...
  fnName <- as.character(substitute(fn, env = environment()))

  summaryTable <- summariseAcrossWindows(qseaSet,
                                         regionsToOverlap = regionsToOverlap,
                                         fn = fn,
                                         suffix = suffix, addSampleTable = FALSE,
                         normMethod = normMethod,
                         naMethod = naMethod,
                         minEnrichment = minEnrichment,
                         fnName = fnName)

  qseaSet <- qseaSet %>%
    left_join(summaryTable, by = "sample_name")

  return(qseaSet)

}


#' Summarise signal by genomic context
#'
#' Annotate windows, then summarise signal and counts by CpG landscape
#' (Island/Shore/Shelf/Open Sea) and short genomic annotation
#' (e.g., Promoter/Exon/Intron/Intergenic).
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing windows and signal.  
#'   **Default:** none (must be supplied).
#'
#' @param cutoff `numeric(1)`.  
#'   Threshold used to call a window “over cutoff” for counting.  
#'   **Default:** `1`.
#'
#' @param normMethod `character(1)`.  
#'   Measure to summarise. One of `"nrpm"` or `"beta"`.  
#'   **Default:** `"nrpm"`.
#'
#' @param minEnrichment `integer(1)`.  
#'   Minimum reads per window required for a non-`NA` beta (relevant when
#'   `normMethod = "beta"`).  
#'   **Default:** `3`.
#'
#' @details
#' Windows are annotated via [annotateWindows()] (which can add CpG landscape
#' and a concise `shortAnno` label). The function then groups by `sample_name`,
#' `landscape`, and `shortAnno` to compute:
#' - `nWindows`: number of windows in the group,  
#' - `sum`: sum of the selected measure across those windows,  
#' - `nOverCutoff`: number of windows with value `>= cutoff`.  
#' When using beta, windows failing `minEnrichment` are treated as `NA` before
#' summarisation.
#' 
#' This function requires a transcript database and OrgDb to annotate windows (via
#' [annotateWindows()]). Either:
#' - set them globally with [setMesaGenome()] (recommended for GRCh38/hg38), or
#' - provide compatible `TxDb`/`annoDb` packages yourself.
#'
#'
#' @return
#' A tibble with columns:
#' `sample_name`, `landscape`, `shortAnno`, `nWindows`, `sum`, `nOverCutoff`,
#' followed by sample-table metadata joined from `qseaSet`.
#'
#' @seealso
#' [annotateWindows()], [summariseAcrossWindows()], [getDataTable()]
#'
#' @family annotation-summaries
#'
#' @examples
#' # Ensure annotation defaults are available (GRCh38/hg38)
#' setMesaGenome("hg38")
#'
#' data(exampleTumourNormal, package = "mesa")
#'
#' # NRPM-based context summary
#' exampleTumourNormal %>%
#'   getGenomicFeatureDistribution(cutoff = 1, normMethod = "nrpm") %>%
#'   head()
#'
#' # Beta-based summary (apply a minimum enrichment for non-NA betas)
#' exampleTumourNormal %>%
#'   getGenomicFeatureDistribution(cutoff = 0.7, normMethod = "beta", minEnrichment = 3) %>%
#'   head()
#'
#' @export
getGenomicFeatureDistribution <- function(qseaSet, cutoff = 1 , normMethod = "nrpm", minEnrichment = 3){
  
  #TODO: This requires the annotateWindows parameters to be exposed to not require setMesaTxDb and setMesaAnnoDb.
  #TODO: Ensure this works when landscape is not present, as that is optional output from annotateWindows
  
  temp <- qseaSet %>%
    getDataTable(normMethod = normMethod,
                 minEnrichment = minEnrichment,
                 addMethodSuffix = TRUE) %>%
    annotateWindows()

  nWindows <- temp %>%
    dplyr::group_by(landscape) %>%
    dplyr::summarise(dplyr::across(tidyselect::matches("nrpm|beta"),~sum(!is.na(.), na.rm = TRUE))) %>%
    tidyr::pivot_longer(tidyselect::matches("nrpm|beta"), names_to = "sample_name", values_to = "nWindows") %>%
    dplyr::mutate(sample_name = stringr::str_remove(sample_name,"_nrpm$|_beta$"))

  sumData <- temp %>%
    dplyr::group_by(landscape) %>%
    dplyr::summarise(dplyr::across(tidyselect::matches("nrpm|beta"),~sum(., na.rm = TRUE))) %>%
    tidyr::pivot_longer(tidyselect::matches("nrpm|beta"), names_to = "sample_name", values_to = "sum") %>%
    dplyr::mutate(sample_name = stringr::str_remove(sample_name,"_nrpm$|_beta$"))

  overCutoffData <- temp %>%
    dplyr::group_by(landscape) %>%
    dplyr::summarise(dplyr::across(tidyselect::matches("nrpm|beta"),~sum(.>= cutoff, na.rm = TRUE))) %>%
    tidyr::pivot_longer(tidyselect::matches("nrpm|beta"), names_to = "sample_name", values_to = "nOverCutoff") %>%
    dplyr::mutate(sample_name = stringr::str_remove(sample_name,"_nrpm$|_beta$"))

  nWindowsShortAnno <- temp %>%
    dplyr::group_by(shortAnno) %>%
    dplyr::summarise(dplyr::across(tidyselect::matches("nrpm|beta"),~sum(!is.na(.), na.rm = TRUE))) %>%
    tidyr::pivot_longer(tidyselect::matches("nrpm|beta"), names_to = "sample_name", values_to = "nWindows") %>%
    dplyr::mutate(sample_name = stringr::str_remove(sample_name,"_nrpm$|_beta$"))

  sumDataShortAnno <- temp %>%
    dplyr::group_by(shortAnno) %>%
    dplyr::summarise(dplyr::across(tidyselect::matches("nrpm|beta"),~sum(.,na.rm = TRUE)), nWindows = dplyr::n()) %>%
    tidyr::pivot_longer(tidyselect::matches("nrpm|beta"), names_to = "sample_name", values_to = "sum") %>%
    dplyr::mutate(sample_name = stringr::str_remove(sample_name,"_nrpm$|_beta$"))

  overCutoffShortAnno <- temp %>%
    dplyr::group_by(shortAnno) %>%
    dplyr::summarise(dplyr::across(tidyselect::matches("nrpm|beta"),~sum(.>= cutoff, na.rm = TRUE))) %>%
    tidyr::pivot_longer(tidyselect::matches("nrpm|beta"), names_to = "sample_name", values_to = "nOverCutoff") %>%
    dplyr::mutate(sample_name = stringr::str_remove(sample_name,"_nrpm$|_beta$"))

  sumData %>%
    dplyr::left_join(overCutoffData) %>%
    dplyr::left_join(nWindows) %>%
    dplyr::full_join(sumDataShortAnno %>%
                       dplyr::left_join(overCutoffShortAnno) %>%
                       dplyr::left_join(nWindowsShortAnno)
    ) %>%
    dplyr::select(sample_name, landscape, shortAnno, nWindows, sum, nOverCutoff) %>%
    dplyr:: left_join(qseaSet %>% qsea::getSampleTable()) %>%
    return()

}

setMethod('getSampleNames', 'data.frame',function(object){stop("getSampleNames is not defined on a data frame, only on a qseaSet.")})


#' Sample names per group (order-preserving)
#'
#' Build a named list mapping each `group` to the vector of `sample_name`s
#' in that group, preserving the current sample order.
#'
#' @param qseaSet `qseaSet`.  
#'   Object providing the sample table (must contain `sample_name` and `group`).  
#'   **Default:** none (must be supplied).
#'
#' @return
#' A named `list` of `character` vectors, where each list name is a group and
#' each element is the ordered vector of sample names in that group.
#'
#' @details
#' Groups are taken from the `group` column of the `sampleTable`. The order of
#' samples within each group reflects their current order in the `qseaSet`.
#' Useful for per-group plotting or aggregation utilities.
#'
#' @seealso
#' [qsea::getSampleTable()], [pull.qseaSet()], [select.qseaSet()]
#'
#' @family table-helpers
#' @keywords internal
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Internal helper (not exported): inspect structure of the mapping
#' exampleTumourNormal %>%
#'   mesa:::getSampleGroups2() 
#'   
getSampleGroups2 <- function(qseaSet){
  qseaSet %>%
    qsea::getSampleTable() %>%
    pull(group) %>%
    unique() %>%
    rlang::set_names(., nm = .) %>%
    purrr::map(function(x){
      qseaSet %>%
        qsea::getSampleTable() %>%
        filter(group == !!x) %>%
        pull(sample_name)
    })
}


#' Extract per-window data (counts / NRPM / beta)
#'
#' Build a per-window table for one or more normalisation methods, returning
#' values per **sample** or **group means**.
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'   **Default:** none (must be supplied).
#'
#' @param normMethod `character()`.  
#'   One or more of `"counts"`, `"nrpm"`, `"beta"`.  
#'   **Default:** `"nrpm"`.
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, average replicates by the `group` column and return group means;
#'   if `FALSE`, return per-sample values.  
#'   **Default:** `FALSE`.
#'
#' @param minEnrichment `integer(1)`.  
#'   Minimum reads per window required for a non-`NA` beta (applies only when
#'   `normMethod` includes `"beta"`).  
#'   **Default:** `3`.
#'
#' @param addMethodSuffix `logical(1)`.  
#'   Keep a method suffix in column names (e.g., `Sample1_nrpm`,
#'   `Group_beta_means`). If multiple methods are requested, suffixes are kept
#'   regardless to avoid collisions.  
#'   **Default:** `FALSE`.
#'
#' @param verbose `logical(1)`.  
#'   Print progress/messages.  
#'   **Default:** `TRUE`.
#'
#' @details
#' The result is a **window × sample (or group)** table. Window coordinates and
#' any available region metadata (e.g., `seqnames`, `start`, `end`, densities)
#' are kept alongside the requested measures. Column naming: when a single
#' method is requested and `addMethodSuffix = FALSE`, columns are bare sample /
#' group names; otherwise a `_{method}` (and for group means `_{method}_means`)
#' suffix is appended.
#'
#' For `"beta"`, windows with reads `< minEnrichment` are set to `NA` before
#' aggregation. Use [removeNormMethodSuffix()] afterwards if you want to strip
#' method tags in wide tables.
#'
#' @return
#' A tibble/data.frame with one row per genomic window and one column per
#' sample (or group), plus window metadata columns.
#'
#' @seealso
#' [getCountTable()], [getNRPMTable()], [getBetaTable()], [removeNormMethodSuffix()]
#'
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Default: NRPM per window (first few columns)
#' exampleTumourNormal %>%
#'   getDataTable() %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Multiple methods together (suffixes retained automatically)
#' exampleTumourNormal %>%
#'   getDataTable(normMethod = c("nrpm","beta"), minEnrichment = 3) %>%
#'   dplyr::select(1:8) %>%
#'   head()
#'
#' # Group means instead of individual samples
#' exampleTumourNormal %>%
#'   getDataTable(useGroupMeans = TRUE, normMethod = "nrpm") %>%
#'   dplyr::select(1:6) %>%
#'   head()
#'
#' # Keep (or later remove) the method suffix
#' exampleTumourNormal %>%
#'   getDataTable(normMethod = "nrpm", addMethodSuffix = TRUE) %>%
#'   names() %>%
#'   head()
#'
#' @export
getDataTable <- function(qseaSet, normMethod = "nrpm", useGroupMeans = FALSE, minEnrichment = 3, addMethodSuffix = FALSE, verbose = TRUE){

  if(!is.qseaSet(qseaSet)){
    stop("Please provide a qseaSet as the first argument.")
  }
  
  if(qseaSet %>% qsea::getRegions() %>% length() == 0){
    stop("Attempting to get data values for a qseaSet with no remaining windows.")
  }
  
  if(useGroupMeans){
    if(verbose){message(glue::glue("Generating table of {normMethod} values for {qseaSet %>% qsea::getRegions() %>% length()} regions across {qseaSet %>% getSampleGroups2() %>% length()} sample groups."))}
    tab <- qseaSet %>%
      qsea::makeTable(groupMeans =  getSampleGroups2(.), 
                      norm_methods = normMethod, 
                      minEnrichment = minEnrichment,
                      verbose = FALSE) %>% #don't use makeTable's messages as we have a different one above.
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble()

    if(length(normMethod) == 1 & !addMethodSuffix) {
      tab <- tab %>%
        dplyr::rename_with(~ stringr::str_replace_all(.x, glue::glue("_{normMethod}_means"), ""))
    } else {
      tab <- tab %>%
        dplyr::rename_with(~ stringr::str_replace_all(.x, "_means$", ""))
    }

    return(tab)

  } else {
    if(verbose){message(glue::glue("Generating table of {normMethod} values for {qseaSet %>% qsea::getRegions() %>% length()} regions across {qseaSet %>% getSampleGroups2() %>% length()} samples."))}
    tab <- qseaSet %>%
      qsea::makeTable(samples = qsea::getSampleNames(.), 
                      norm_methods = normMethod, 
                      minEnrichment = minEnrichment,
                      verbose = FALSE) %>% #don't use makeTable's messages as we have a different one above.
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble()

    if(length(normMethod) == 1 & !addMethodSuffix) {
      tab <- tab %>%
        dplyr::rename_with(~ stringr::str_replace_all(.x, glue::glue("_{normMethod}"), ""))
    }

    return(tab)
  }

}


#' Write bigWig tracks per sample or group
#'
#' Export bigWig files with per-window scores for each sample (or group means).
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'   **Default:** none (must be supplied).
#'
#' @param folderName `character(1)`.  
#'   Output directory where `.bw` files will be written (created if missing).  
#'   **Default:** none (must be supplied).
#'
#' @param normMethod `character(1)`.  
#'   Normalisation/measure to export. Common choices: `"nrpm"` or `"beta"`.  
#'   **Default:** `"nrpm"`.
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, export group-mean tracks instead of per-sample tracks.  
#'   **Default:** `FALSE`.
#'
#' @param naVal `numeric(1)`.  
#'   Value used to replace `NA` scores prior to export (useful for beta).  
#'   **Default:** `0`.
#'
#' @details
#' For each sample (or group), the function pairs the per-window scores from
#' `qseaSet` with the corresponding window genomic coordinates and writes a
#' bigWig via **rtracklayer**. When exporting beta, windows failing the minimum
#' read threshold used to compute beta may be `NA`; these are replaced by
#' `naVal` before writing.
#'
#' **File naming.** Output files are typically named using the sample (or group)
#' name and the method, e.g. `"<sample>_<method>.bw"` or
#' `"<group>_<method>_means.bw"` when `useGroupMeans = TRUE`.
#'
#' **Seqinfo.** bigWig export requires chromosome lengths; ensure the window
#' GRanges in `qseaSet` carries valid `seqinfo` (this is the case for the
#' example data on GRCh38).
#'
#' @return
#' (Invisibly) `NULL`. Files are written to `folderName`.
#'
#' @seealso
#' [getDataTable()]
#'
#' @family io
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Per-sample NRPM bigWigs
#' td <- tempfile()
#' dir.create(td)
#' exampleTumourNormal %>%
#'   writeBigWigs(folderName = td, normMethod = "nrpm", useGroupMeans = FALSE)
#' list.files(td, pattern = "\\\\.bw$")
#'
#' # Group-mean beta bigWigs (replace NA beta with 0)
#' td2 <- tempfile()
#' dir.create(td2)
#' exampleTumourNormal %>%
#'   writeBigWigs(folderName = td2, normMethod = "beta", useGroupMeans = TRUE, naVal = 0)
#' list.files(td, pattern = "\\\\.bw$")
#'
#' @export
writeBigWigs <- function(qseaSet, folderName, normMethod = "nrpm", useGroupMeans = FALSE, naVal = -1){

  dir.create(folderName, showWarnings = FALSE)

  dataTable <- qseaSet %>%
    getDataTable(normMethod = normMethod, useGroupMeans = useGroupMeans)

  if(!useGroupMeans) {
    mapNames <- qseaSet %>%
      qsea::getSampleNames()

  } else {
    mapNames <- qseaSet %>%
      getSampleGroups2()

  }

  mapNames %>%
    purrr::walk(function(x){

      message(glue::glue("Writing bigWig track to {folderName}/{x}_{normMethod}.bw"))

      gr <- dataTable %>%
        dplyr::rename(score = !!x) %>%
        dplyr::select(seqnames, start, end, score) %>%
        dplyr::mutate(score = tidyr::replace_na(score, naVal)) %>%
        plyranges::as_granges()

      GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(qseaSet %>% qsea::getRegions())

      gr %>%
        plyranges::write_bigwig(glue::glue("{folderName}/{x}_{normMethod}.bw"))

    })
}


#' Set library factors to 1
#'
#' Overwrite all library-factor values in a `qseaSet` (and its `sampleTable`)
#' with `1`. Useful when you want to avoid library scaling (e.g., TMM) and
#' keep beta unchanged; NRPM may change depending on downstream use.
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'   **Default:** none (must be supplied).
#'
#' @details
#' Setting library factors to `1` has **no effect on beta** (fractional measure)
#' but can affect **NRPM** since scaling factors are neutralised. This is
#' akin to calling `qsea::addLibraryFactors(qseaSet, factors = 1)` and
#' synchronising the `sampleTable`.
#'
#' @return
#' A `qseaSet` with all library factors set to `1` and `sampleTable` updated.
#'
#' @seealso
#' [addLibraryInformation()], [getNRPMTable()], [qsea::addLibraryFactors()]
#'
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Show current factors then reset to 1 and show again
#' exampleTumourNormal %>%
#'   addLibraryInformation() %>%
#'   pull(library_factor) %>%
#'   head()
#'
#' @export
removeLibraryFactors <- function(qseaSet){
  qseaSet <- qseaSet %>%
    qsea::addLibraryFactors(1)

  if("library_factor" %in% colnames(qsea::getSampleTable(qseaSet))){
    qseaSet <- qseaSet %>%
      dplyr::mutate(library_factor = 1)
  }
return(qseaSet)
}
