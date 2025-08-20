# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("Welcome to my methtools package!")
# }


# Internal: memoise on load
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  getCGPositions <<- memoise::memoise(getCGPositions)
}


#' Check if an object is a qseaSet
#'
#' This function checks that an object is a qseaSet.
#'
#' @param x The object
#' 
#' @return Logical scalar: `TRUE` if `x` is a `qseaSet`, `FALSE` otherwise.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("qsea", quietly = TRUE)) {
#'   if (system.file("data", "exampleTumourNormal.rda", package = "mesa") != "") {
#'     data("exampleTumourNormal", package = "mesa")
#'     is.qseaSet(exampleTumourNormal)  # TRUE
#'   }
#'   is.qseaSet(iris)                    # FALSE
#' }
#' }
#' @export
is.qseaSet <- function(x){
  return(inherits(x,"qseaSet"))
  }


#' Set or get an Ensembl/Biomart (or other) handle on a qseaSet
#'
#' These functions store and retrieve a “mart” object (or string handle) inside a `qseaSet`
#' for downstream annotation tasks.
#'
#' @name setMart
#' @param object A `qseaSet`.
#' @param mart   An object or string identifying the mart to use (e.g., an Ensembl/BioMart
#'   connection or a label you interpret elsewhere).
#'
#' @return
#' - `setMart()` returns the updated `qseaSet` (with `@parameters$mart` set).
#' - `getMart()` returns the stored mart value.
#'
#' @seealso [annotateWindows()], [setMesaGenome()], [setMesaTxDb()], [setMesaAnnoDb()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("qsea", quietly = TRUE)) {
#'   if (system.file("data", "exampleTumourNormal.rda", package = "mesa") != "") {
#'     data("exampleTumourNormal", package = "mesa")
#'     qs <- setMart(exampleTumourNormal, mart = "ENSEMBL_110")
#'     getMart(qs)
#'   }
#' }
#' }
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
#' @details
#' - **Inputs:** `dataTable` may be a `GRanges`, or a data frame coercible by
#'   `qseaTableToChrGRanges()` (must include `seqnames`, `start`, `end`).
#' - **Genome defaults (GRCh38):** When `genome` is `"hg38"`/`"GRCh38"` and `TxDb`
#'   or `annoDb` are `NULL`, the function uses
#'   `TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene`
#'   and `"org.Hs.eg.db"` if installed. For `CpGislandsGR` and `FantomRegionsGR`,
#'   mesa defaults (`hg38CpGIslands`, `FantomRegions`) are used if available.
#' - **Seqname style:** The output removes the `"chr"` prefix (`seqnames` become `"1"…"22","X","Y"`).
#' - **Output columns:** Adds `shortAnno` (annotation without the trailing parenthetical),
#'   optional CpG landscape (`landscape` with Island/Shore/Shelf/Open Sea) and
#'   `inFantom` (overlap count with FANTOM enhancers). Width and strand are dropped.
#'
#' @param dataTable A data frame coercible to `GRanges`, or a `GRanges` directly.
#' @param genome Genome string guiding defaults (currently `"hg38"`/`"GRCh38"` supported).
#'   If `NULL`, a fully specified `TxDb`/`annoDb` must be supplied.
#' @param TxDb   A TxDb object (unquoted) or a string of the form
#'   `"PkgName::ObjectName"`. If `NULL` and `genome` is GRCh38, a default TxDb is used.
#' @param annoDb A string naming a Bioconductor OrgDb package (e.g., `"org.Hs.eg.db"`).
#'   If `NULL` and `genome` is GRCh38, `"org.Hs.eg.db"` is used if installed.
#' @param CpGislandsGR Optional `GRanges` of CpG islands (if `NULL` and GRCh38, uses `mesa::hg38CpGIslands` if available).
#' @param FantomRegionsGR Optional `GRanges` of FANTOM enhancers (if `NULL` and GRCh38, uses `mesa::FantomRegions` if available).
#'
#' @return A tibble with the input windows augmented by ChIPseeker annotations and,
#' if provided/available, CpG island landscape and FANTOM overlap counts.
#'
#' @seealso
#' [ChIPseeker::annotatePeak()], [setMesaGenome()], [setMesaTxDb()], [setMesaAnnoDb()],
#' [liftOverHg19()]
#'
#' @examples
#' \donttest{
#' ok <- requireNamespace("ChIPseeker", quietly = TRUE) &&
#'       requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
#'       requireNamespace("org.Hs.eg.db", quietly = TRUE)
#' if (ok) {
#'   # Minimal toy GRanges
#'   gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1e6 + 0:4 * 1000, width = 200))
#'
#'   # Use GRCh38 defaults (TxDb/org.Hs.eg.db; mesa CpG/FANTOM if available)
#'   setMesaGenome("hg38")
#'   ann <- annotateWindows(gr)
#'   head(ann)
#'
#'   # Data-frame input also works (coerced internally)
#'   df <- data.frame(seqnames = "chr1", start = 2e6, end = 2e6 + 199)
#'   annotateWindows(df)
#' }
#' }
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
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }

  if(genome  %in% c("hg38","GRCh38") && is.null(annoDb)) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop(
        "Package \"org.Hs.eg.db\" must be installed to use this function. Please install and run again.",
        call. = FALSE
      )
    }
    annoDb = "org.Hs.eg.db"
    }

  if(is.null(annoDb) && is.null(genome)) {
    stop("Please specify a annoDb or genome, this can be set globally using setMesaannoDb and/or setMesaGenome")
  }
  
  if(genome  %in% c("hg38","GRCh38") && is.null(CpGislandsGR)) { CpGislandsGR = mesa::hg38CpGIslands }

  if(genome  %in% c("hg38","GRCh38") && is.null(FantomRegionsGR)) { FantomRegionsGR = mesa::FantomRegions %>% plyranges::as_granges()}

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
#' Filter genomic windows based on a summary function applied to selected
#' samples. For example, retain windows with median normalized counts above
#' a threshold, or where the minimum beta value is below 0.5.
#'
#' @param qseaSet A `qseaSet` object.
#' @param fn A function applied row-wise across selected samples (e.g. `median`,
#'   `min`, `max`).
#' @param threshold `numeric(1)` Threshold value applied to `fn` results.
#' @param aboveThreshold `logical(1)` If `TRUE`, retain windows with values
#'   `>= threshold`; otherwise, retain windows `< threshold`.
#' @param samples `character()` Sample names to consider, or a string pattern
#'   matched against sample names. Defaults to all samples.
#' @param normMethod `character(1)` Normalization method, one of `"nrpm"`,
#'   `"beta"`, or `"counts"`.
#' @param useGroupMeans `logical(1)` If `TRUE`, aggregate replicate samples
#'   by their `group` in the sample table.
#'
#' @return A filtered `qseaSet` containing only the selected windows.
#'
#' @seealso [getDataTable()], [filterByOverlaps()]
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#'
#' ## Keep windows with median nrpm > 1 across tumour samples
#' subsetWindowsBySignal(qs, fn = median, threshold = 1,
#'                       aboveThreshold = TRUE, samples = "LUAD")
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

  #TODO think about catching the warnings more specifically. Mostly we want to catch the "no non-missing arguments to max; returning -Inf"
  dataTable <- suppressWarnings(dataTable %>% dplyr::mutate(fnValue = apply(dplyr::pick(tidyselect::all_of(samples)), 1, fn, na.rm = TRUE)))

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
#' Identify windows with significantly more reads than expected under a
#' Poisson background model, given the total reads per sample and number of
#' windows in the genome.
#'
#' @param qseaSet A `qseaSet` object.
#' @param keepAbove `logical(1)` If `TRUE`, retain windows above background;
#'   if `FALSE`, remove them.
#' @param samples `character()` Vector of sample names to test, or a string
#'   pattern matched against sample names. Defaults to all samples.
#' @param numWindows `integer(1)` Number of windows in the genome. If `NULL`,
#'   estimated automatically.
#' @param recalculateNumWindows `logical(1)` Whether to recalculate the number
#'   of windows via [qsea::createQseaSet()].
#' @param fdrThres `numeric(1)` FDR threshold for significance.
#' @param numAbove `integer(1)` Minimum number of samples that must exceed
#'   background in a window to retain/drop it.
#'
#' @return A filtered `qseaSet` containing only windows passing the criteria.
#'
#' @seealso [subsetWindowsBySignal()]
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#'
#' ## Keep only windows above Poisson background in at least 2 tumour samples
#' subsetWindowsOverBackground(qs, keepAbove = TRUE,
#'                             samples = "LUAD", numAbove = 2)
#'
#' @importFrom rlang :=
#' @export
subsetWindowsOverBackground <- function(qseaSet, keepAbove = FALSE, samples = NULL, numWindows = NULL,
                                        recalculateNumWindows = TRUE, fdrThres = 0.01, numAbove = 1){


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

    if (recalculateNumWindows) {

      numWindows <- suppressWarnings(qsea::createQseaSet(sampleTable = qsea::getSampleTable(qseaSet),
                                                   BSgenome = qseaSet@parameters$BSgenome,
                                                   chr.select = qsea::getChrNames(qseaSet),
                                                   window_size = qsea::getWindowSize(qseaSet))
      ) %>%
        qsea::getRegions() %>%
        length()
    } else {
      numWindows <- nrow(countMat)
    }
  }

  fdrMat <- purrr::map_dfc(samples,
                    function(x){

                      totalNumReads <- qseaSet@libraries$file_name[x, "total_fragments"]
                      lambda <- totalNumReads/numWindows
                      pvals <- stats::ppois(countMat[,x] - 1, lambda, lower.tail = FALSE)
                      fdrvals <- stats::p.adjust(pvals, method = "fdr") %>%
                        tibble::enframe(name = "window") %>%
                        dplyr::rename(!!x := value) %>%
                        dplyr::select(-window)

                    }
  )

  if (keepAbove) {
    windowsToKeep <- fdrMat %>%
      {. <= fdrThres} %>%
      rowSums() %>%
      {. >= numAbove} %>%
      which()
  } else {
    windowsToKeep <- fdrMat %>%
      {. <= fdrThres} %>%
      rowSums() %>%
      {. < numAbove} %>%
      which()
  }

  message(glue::glue("Removing {nrow(fdrMat) - length(windowsToKeep)} windows based on {length(samples)} samples, {length(windowsToKeep)} remaining"))

  return(invisible(filterByOverlaps(qseaSet, windowsToKeep)))

}


#' Downsample reads in a qseaSet
#'
#' Randomly subsample reads to a fixed depth across all samples in a `qseaSet`.
#' This is useful for equalizing library sizes prior to comparison.
#'
#' @param qseaSet A `qseaSet` object.
#' @param nReads `integer(1)` Target number of reads to retain per sample.
#'
#' @return A `qseaSet` with downsampled counts. Library metadata
#' (`valid_fragments`, `offset`, `library_factor`) are updated accordingly.
#'
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#'
#' ## Downsample to 1e5 reads per sample
#' ds <- downSample(qs, nReads = 1e5)
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

        sample(rep(1:length(vec), vec), replace = FALSE, size = nReads) %>%
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


#' Convert qsea beta values to array-like format
#'
#' Extract beta values for probes from a methylation array, suitable for use
#' with tools expecting array-style beta matrices. Currently supports
#' Illumina Infinium450k probes on GRCh38.
#'
#' @param qseaSet A `qseaSet` object.
#' @param arrayDetails Either a recognized string (currently `"Infinium450k"`)
#'   or a `GRanges` object with an `ID` column of probe identifiers.
#'
#' @return A `data.frame` with rows as probes (including `ID` column) and
#'   columns as samples with corresponding beta values.
#'
#' @seealso [getDataTable()]
#' @family window-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' ## Using built-in probe set
#' beta_tab <- convertToArrayBetaTable(exampleTumourNormal,
#'                                     arrayDetails = "Infinium450k")
#'
#' ## Using a GRanges probe object
#' beta_tab2 <- convertToArrayBetaTable(exampleTumourNormal,
#'                                      arrayDetails = mesa::hg38_450kArrayGR)
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
#' For each sample in a `qseaSet`, compute the proportion of **windows** with
#' counts ≥ `numCountsNeeded` that overlap `windowsToConsider`, relative to all
#' windows with counts ≥ `numCountsNeeded` in that sample.
#'
#' @details
#' This function operates on the window-by-sample count matrix from
#' `qsea::getCounts(qseaSet)`. It first converts counts to a logical indicator
#' (count ≥ `numCountsNeeded`) and sums these indicators per sample:
#' - `initialOverBackNum`: number of windows meeting the threshold genome-wide.
#' - `afterOverBackNum`: number of windows meeting the threshold **within**
#'   `windowsToConsider` (via `filterByOverlaps()`).
#' The reported `fraction` is `afterOverBackNum / initialOverBackNum`.
#'
#' Note: despite the function name, this is a **window-based** fraction using a
#' count threshold, not a direct fraction of raw reads. Use accordingly.
#'
#' @param qseaSet A `qseaSet` object.
#' @param windowsToConsider A `GRanges` of regions to consider/overlap.
#' @param numCountsNeeded Integer; minimum reads per window for it to be counted.
#'
#' @return
#' A tibble with one row per sample containing:
#' - `sample_name`
#' - `initialOverBackNum` — windows ≥ threshold genome-wide
#' - `afterOverBackNum` — windows ≥ threshold within `windowsToConsider`
#' - `fraction` — `afterOverBackNum / initialOverBackNum`
#' followed by columns from the sample table (via a left join on `sample_name`).
#'
#' @seealso
#' [qsea::getCounts()], [qsea::getSampleTable()], [filterByOverlaps()], [hg38UltraStableProbes]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("qsea", quietly = TRUE)) {
#'   # Example qseaSet shipped with mesa (if available)
#'   if (system.file("data", "exampleTumourNormal.rda", package = "mesa") != "") {
#'     data("exampleTumourNormal", package = "mesa")
#'
#'     # Use shipped ultra-stable regions (a GRanges) as windowsToConsider
#'     if (system.file("data", "hg38UltraStableProbes.rda", package = "mesa") != "") {
#'       data("hg38UltraStableProbes", package = "mesa")
#'       res <- calculateFractionReadsInGRanges(
#'         exampleTumourNormal,
#'         hg38UltraStableProbes,
#'         numCountsNeeded = 5
#'       )
#'       head(res[, c("sample_name", "fraction")])
#'     }
#'   }
#' }
#' }
#' @export
calculateFractionReadsInGRanges <- function(qseaSet, windowsToConsider, numCountsNeeded) {
  initialReadTotals <- qseaSet %>%
    qsea::getCounts() %>%
    {. >= numCountsNeeded } %>%
    colSums()

  afterSubsetReadTotals <- qseaSet %>%
    filterByOverlaps(windowsToConsider) %>%
    qsea::getCounts() %>%
    {. >= numCountsNeeded } %>%
    colSums()

  tibble::tibble(sample_name = names(initialReadTotals),
         initialOverBackNum = initialReadTotals,
         afterOverBackNum = afterSubsetReadTotals,
         fraction = afterOverBackNum/initialOverBackNum) %>%
    dplyr::left_join(qsea::getSampleTable(addLibraryInformation(qseaSet))) %>%
    return()

}


#' Create a heatmap annotation data.frame (individual samples)
#'
#' Generate an annotation data.frame for use with heatmaps, 
#' containing metadata columns from the `qseaSet` sample table. 
#' Each row corresponds to an individual sample.
#'
#' @param qseaSet A [qsea::qseaSet] object.
#' @param ... One or more unquoted column names from the sample table 
#'   to include in the annotation data.frame.
#'
#' @return A `data.frame` with rownames set to `sample_name` and 
#'   the requested metadata columns.
#'
#' @seealso [getAnnotationDataFrame()] for group-level annotations.
#' @family heatmap-annotation
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' ann <- getAnnotationDataFrameIndividual(qs, group)
#' head(ann)
#' 
getAnnotationDataFrameIndividual <- function(qseaSet, ...){

  #TODO Can this function be removed?
  if (!("valid_fragments" %in% colnames(qsea::getSampleTable(qseaSet)))) {
    qseaSet <- qseaSet %>%
      addLibraryInformation()
  }

  qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::distinct(sample_name,...) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("sample_name") %>%
    return()
}


#' Create a heatmap annotation data.frame (groups)
#'
#' Generate an annotation data.frame for use with heatmaps, 
#' aggregating sample information at the group level.
#'
#' @param qseaSet A [qsea::qseaSet] object.
#' @param ... One or more unquoted column names from the sample table 
#'   to include in the annotation data.frame.
#'
#' @return A `data.frame` with rownames set to `group` and 
#'   the requested metadata columns.
#'
#' @seealso [getAnnotationDataFrameIndividual()] for sample-level annotations.
#' @family heatmap-annotation
#' 
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' ann_grp <- getAnnotationDataFrame(qs, disease)
#' head(ann_grp)
#' 
getAnnotationDataFrame <- function(qseaSet, ...){

  #TODO Can this function be removed?
  if (!("valid_fragments" %in% colnames(qsea::getSampleTable(qseaSet)))) {
    qseaSet <- qseaSet %>%
      addLibraryInformation()
  }

  qseaSet %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(total_fragments = mean(total_fragments),
           relH = mean(relH)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(group,...) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("group") %>%
    return()
}


#' Remove normalisation suffix from column names
#'
#' Drop the `_{normMethod}` (and optional `_means`) suffix from sample or
#' group mean columns in a wide data table (e.g., from [getDataTable()]).
#'
#' @param dataTable data.frame with window x sample/group values.
#' @param normMethod character(1) normalisation method (e.g., `"beta"`, `"nrpm"`).
#'
#' @return data.frame with cleaned column names.
#' @family table-helpers
#'
#' @examples
#' df <- data.frame(A_beta = 1:3, B_beta_means = 4:6)
#' removeNormMethodSuffix(df, "beta")
#' 
#' @export
removeNormMethodSuffix <- function(dataTable, normMethod) {
  dplyr::rename_with(dataTable, ~ stringr::str_remove(.x, glue::glue("_{normMethod}(_means)?$")))
}


#' Count windows above a cutoff
#'
#' For a set of windows, count how many exceed a cutoff per sample.
#'
#' @param qseaSet qsea::qseaSet
#' @param GRanges GenomicRanges::GRanges windows to filter on.
#' @param samples character() sample names to include (default: all).
#' @param normMethod character(1) normalisation method (`"nrpm"` or `"beta"`).
#' @param cutoff numeric(1) threshold to count windows above.
#'
#' @return tibble with columns: `sample_name`, `numOverCutoff`, `totalWindowsUsed`,
#'   plus sample table columns joined from the `qseaSet`.
#' @seealso [getDataTable()], [summariseAcrossWindows()]
#' @family window-summaries
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' gr <- qsea::getRegions(qs)[1:100]
#' out <- countWindowsAboveCutoff(qs, gr, cutoff = 1, normMethod = "nrpm")
#' head(out)
#' 
#' @export
countWindowsAboveCutoff <- function(qseaSet, GRanges, samples = NULL,
                                   cutoff = 0, normMethod = "nrpm"){

  if (is.null(samples)) {
    samples <- qsea::getSampleNames(qseaSet)
  }

  clip <- function(x, a, b) {
    a + (x - a > 0) * (x - a) - (x - b > 0) * (x - b)
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
#' Create a wide table suitable for ML/visualisation where rows are samples and
#' columns are windows, using a single normalisation method.
#'
#' @param qseaSet qsea::qseaSet
#' @param normMethod character(1) normalisation method (e.g., `"nrpm"`).
#' @param ... Optional columns from the sample table to append.
#'
#' @return tibble with one row per sample and one column per genomic window.
#' @seealso [getDataTable()]
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' tbl <- makeTransposedTable(qs, normMethod = "nrpm", group)
#' dim(tbl)
#' 
#' @export
makeTransposedTable <- function(qseaSet, normMethod = "nrpm", ...){

  qseaSet %>%
    qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = normMethod) %>%
    dplyr::rename_with(~ stringr::str_replace_all(.x, "_nrpm$|_beta$|_counts$", ""))  %>%
    dplyr::select(-CpG_density) %>%
    tidyr::pivot_longer(-c(chr, window_start,window_end), names_to = "sample_name", values_to = "value") %>%
    dplyr::mutate(chr = paste0("chr",chr)) %>%
    tidyr::unite(col = "window", chr, window_start, window_end) %>%
    tidyr::pivot_wider(names_from = window, values_from = value) %>%
    dplyr::left_join(qsea::getSampleTable(qseaSet) %>% dplyr::select(sample_name, ...)) %>%
    dplyr::relocate(tidyselect::matches("^chr"), .after = tidyselect::last_col())

}


#' Get counts per window
#'
#' Convenience wrapper around [getDataTable()] to extract *counts*.
#'
#' @param qseaSet qsea::qseaSet object.
#' @param useGroupMeans logical(1). If `TRUE`, average replicates by the `group`
#'   column instead of returning per-sample values.
#' @param addMethodSuffix logical(1). If `TRUE`, keep the method suffix in
#'   column names (e.g., `Sample1_counts`). The suffix is always kept if
#'   multiple normalisation methods are requested.
#' @param verbose logical(1). Print progress/messages.
#'
#' @return A tibble/data.frame of *counts* (one row per window).
#' @seealso [getDataTable()]
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' head(getCountTable(qs))
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
#' Convenience wrapper around [getDataTable()] to extract *NRPM* values.
#'
#' @param qseaSet qsea::qseaSet object.
#' @param useGroupMeans logical(1). If `TRUE`, average replicates by the `group`
#'   column instead of returning per-sample values.
#' @param addMethodSuffix logical(1). If `TRUE`, keep the method suffix in
#'   column names (e.g., `Sample1_nrpm`). The suffix is always kept if
#'   multiple normalisation methods are requested.
#' @param verbose logical(1). Print progress/messages.
#'
#' @return A tibble/data.frame of *NRPM* values (one row per window).
#' @seealso [getDataTable()]
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' head(getNRPMTable(qs))
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
#' Convenience wrapper around [getDataTable()] to extract *beta* values.
#'
#' @param qseaSet qsea::qseaSet object.
#' @param useGroupMeans logical(1). If `TRUE`, average replicates by the `group`
#'   column instead of returning per-sample values.
#' @param minEnrichment integer(1). Minimum number of reads required for a
#'   window to be considered fully methylated (below this, qsea sets beta to NA).
#' @param addMethodSuffix logical(1). If `TRUE`, keep the method suffix in
#'   column names (e.g., `Sample1_beta`). The suffix is always kept if
#'   multiple normalisation methods are requested.
#' @param verbose logical(1). Print progress/messages.
#'
#' @return A tibble/data.frame of *beta* values (one row per window).
#' @seealso [getDataTable()]
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' head(getBetaTable(qs, minEnrichment = 3))
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
#' Apply a summary function (e.g. mean, median, sd) over windows per sample,
#' optionally restricted to a set of regions, and for one or more methods.
#'
#' @param qseaSet qsea::qseaSet
#' @param regionsToOverlap GRanges or data.frame coercible to GRanges; if `NULL`,
#'   use all regions in the qseaSet.
#' @param fn function summary function (e.g., `mean`).
#' @param suffix character(1) suffix for output columns (e.g., a region label).
#' @param addSampleTable logical whether to join sample table columns.
#' @param normMethod character() one or more of `"nrpm"`, `"beta"`.
#' @param naMethod character(1) `"na.rm"` (default) or `"drop"` rows with NA.
#' @param minEnrichment integer(1) minimum reads for non-NA beta.
#' @param fnName character(1) name of `fn` (auto-detected if `NULL`).
#'
#' @return tibble with one row per sample and summary columns per method.
#' @seealso [addSummaryAcrossWindows()], [getDataTable()]
#' @family window-summaries
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' gr <- qsea::getRegions(qs)[1:500]
#' sm <- summariseAcrossWindows(qs, regionsToOverlap = gr, fn = mean,
#'                              normMethod = c("nrpm","beta"), suffix = "first500")
#' head(sm)
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
    fnName = as.character(substitute(fn, env = environment()))
  }

    #if suffix doesn't start with "_" then add that to the string
    suffix = ifelse(stringr::str_detect(suffix, "^_") | nchar(suffix) == 0 , suffix, paste0("_",suffix))

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
#' @inheritParams summariseAcrossWindows
#'
#' @return qsea::qseaSet with new summary columns in the sample table.
#' @seealso [summariseAcrossWindows()]
#' @family window-summaries
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' gr <- qsea::getRegions(qs)[1:300]
#' qs2 <- addSummaryAcrossWindows(qs, regionsToOverlap = gr, fn = mean,
#'                                normMethod = "nrpm", suffix = "first300")
#' head(qsea::getSampleTable(qs2))
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
  fnName = as.character(substitute(fn, env = environment()))

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
#' (island/shore/shelf) and short genomic annotation (e.g. promoter/exon/intron).
#'
#' @param qseaSet qsea::qseaSet
#' @param cutoff numeric(1) threshold to call window “over cutoff”.
#' @param normMethod character(1) `"nrpm"` or `"beta"`.
#' @param minEnrichment integer(1) minimum reads for non-NA beta.
#'
#' @return tibble with `sample_name`, `landscape`, `shortAnno`,
#'   `nWindows`, `sum`, `nOverCutoff`, and sample metadata columns.
#' @seealso [annotateWindows()], [summariseAcrossWindows()]
#' @family annotation-summaries
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' gfd <- getGenomicFeatureDistribution(qs, cutoff = 1, normMethod = "nrpm")
#' head(gfd)
#' 
#' @export
getGenomicFeatureDistribution <- function(qseaSet, cutoff = 1 , normMethod = "nrpm", minEnrichment = 3){
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
#' Construct a named list mapping each `group` to the vector of `sample_name`s
#' in that group, preserving original ordering.
#'
#' @param qseaSet qsea::qseaSet
#' @return named list of character vectors.
#' @keywords internal
#' @family table-helpers
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' # getSampleGroups2(qs)
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


#' Extract per-window data (counts/NRPM/beta)
#'
#' Generate a per-window table for one or more normalisation methods,
#' for individual samples or group means.
#'
#' @param qseaSet qsea::qseaSet
#' @param normMethod character() one or more of `"counts"`, `"nrpm"`, `"beta"`.
#' @param useGroupMeans logical use group means instead of individual samples.
#' @param minEnrichment integer(1) minimum reads for non-NA beta.
#' @param addMethodSuffix logical keep method suffix in column names.
#' @param verbose logical print messages.
#'
#' @return tibble with one row per window and columns per sample or group.
#' @seealso [getCountTable()], [getNRPMTable()], [getBetaTable()]
#' @family table-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' head(getDataTable(qs, normMethod = "nrpm"))
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
#' @param qseaSet qsea::qseaSet
#' @param folderName character(1) output folder.
#' @param normMethod character(1) normalisation method (e.g., `"nrpm"` or `"beta"`).
#' @param useGroupMeans logical export group means instead of samples.
#' @param naVal numeric(1) value used to replace `NA` scores (for beta).
#'
#' @return (invisible) `NULL`. Files are written to `folderName`.
#' @seealso [getDataTable()]
#' @family io
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' td <- tempfile(); dir.create(td)
#' writeBigWigs(qs, td, normMethod = "nrpm", useGroupMeans = FALSE)
#' list.files(td, pattern = "\\.bw$")
#' }
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
#' Set all library-factor columns in the qseaSet (and sample table) to one.
#'
#' @param qseaSet qsea::qseaSet
#'
#' @return qsea::qseaSet with library factors set to 1.
#' @family table-helpers
#' @export
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' qs2 <- removeLibraryFactors(qs)
#' # qsea::getSampleTable(qs2)$library_factor
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

