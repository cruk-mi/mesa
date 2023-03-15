# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("Welcome to my methtools package!")
# }

.onLoad <- function(libname, pkgname) {
  getCGPositions <<- memoise::memoise(getCGPositions)
}

#' Check if an object is a qseaSet
#'
#' This function checks that an object is a qseaSet.
#'
#' @param x The object
is.qseaSet <- function(x){
  return(inherits(x,"qseaSet"))
  }


setGeneric('setMart', function(object,...) standardGeneric('setMart'))
setMethod('setMart', 'qseaSet', function(object, martName){
  object@parameters$mart = martName
  object
})

setGeneric('getMart', function(object,...) standardGeneric('getMart'))
setMethod('getMart', 'qseaSet', function(object)
  object@parameters$mart
)

#' Add annotation onto a data table
#'
#' This function uses the ChIPseeker::annotatePeak function to determine the closest region to each genomic window provided.
#' Defaults to being hg38, unless
#'
#' @param dataTable A data frame which can be coerced into a GRanges object or a GRanges object directly.
#' @param genome A genome string to set the rest of the parameters (currently only hg38/GRCh38 supported)
#' @param TxDb A TxDb database object (unquoted) to pass to ChIPseeker::annotatePeak
#' @param annoDb A string giving a Bioconductor annotation package, such as "org.Hs.eg.db"
#' @param CpGislandsGR A GRanges object giving locations of CpG islands
#' @param FantomRegionsGR A GRanges object giving Fantom enhancer regions
#' @return A tibble with the data, augmented with ChIPseeker region location and CpG island information.
#' @export
annotateWindows <- function(dataTable, genome = "hg38", TxDb = NULL, annoDb = NULL, CpGislandsGR = NULL,
                         FantomRegionsGR = NULL) {

  if(genome %in% c("hg38","GRCh38") & is.null(TxDb)) {

    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      stop(
        "Package \"TxDb.Hsapiens.UCSC.hg38.knownGene\" must be installed to use this function. Please install and run again.",
        call. = FALSE
      )
    }
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }

  if(genome  %in% c("hg38","GRCh38") & is.null(annoDb)) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop(
        "Package \"org.Hs.eg.db\" must be installed to use this function. Please install and run again.",
        call. = FALSE
      )
    }
    annoDb = "org.Hs.eg.db"
    }

  if(genome  %in% c("hg38","GRCh38") & is.null(CpGislandsGR)) { CpGislandsGR = mesa::hg38CpGIslands }

  if(genome  %in% c("hg38","GRCh38") & is.null(FantomRegionsGR)) { FantomRegionsGR = mesa::FantomRegions %>% plyranges::as_granges()}


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

#' Filter windows from a qseaSet based on expression in a subset of samples
#'
#' This function is designed to take a qseaSet and filter the windows by applying a function to each window, followed by thresholding either windows above or below that cutoff.
#' For instance, we can keep only regions that have a median normalised reads per million value above 1, or those where the minimum beta value is below 0.5.
#' Which samples to apply the filter based on can be specified with the `samples` argument, either as a list of samples or a string.
#'
#' @param qseaSet The qseaSet object.
#' @param fn A function that takes a row of data and returns a single value.
#' @param threshold The value to threshold the `fn` output value on.
#' @param aboveThreshold A logical value indicating whether to keep the windows with `fn` output value above or equal to `threshold` (TRUE), or below `threshold` (FALSE).
#' @param samples A set of sample names to filter on, or a string to match in the sample names.
#' @param normMethod The type of normalisation method to use (e.g. nrpm, beta, counts).
#' @param useGroupMeans Whether to use the group argument of the sample table to collapse replicates.
#' @return A qseaSet object, with only a subset of the windows.
#' @export
subsetWindowsBySignal <- function(qseaSet, fn, threshold, aboveThreshold, samples = NULL, normMethod = "nrpm", useGroupMeans = FALSE){

  fnName <- as.character(substitute(fn, env = environment()))

  if (!useGroupMeans) {
    qseaSamples <- qsea::getSampleNames(qseaSet)
    groupString <- ""
  } else {
    qseaSamples <- names(qsea::getSampleGroups(qseaSet))
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

#' This function takes a qseaSet and finds which windows have more reads than would be expected if the reads followed a Poisson distribution over the whole genome.
#'
#'#' Moved to internal function as of January 2023, think whether want to keep it or not.
#'
#' @param qseaSet The qseaSet object.
#' @param keepAbove Boolean, if true, then keep only the windows that are above the background, if false then remove them.
#' @param samples A vector of samples to remove, or a character string to match and remove.
#' @param numWindows The number of windows in the whole genome.
#' @param recalculateNumWindows Whether to recalculate the number of windows based on the whole genome.
#' @param fdrThres What FDR p-value to filter windows out based on.
#' @param numAbove How many samples need to be above the background level in a window to keep/remove it.
#' @return A qseaSet object, with only a subset of the windows.
#' @importFrom rlang :=
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

#' This function takes a qseaSet and downsamples the reads.
#' @param qseaSet The qseaSet object.
#' @param nReads How many reads to downsample to.
#' @param renormalise Whether to renormalise the resulting qseaSet.
#' @return A qseaSet object with the reads downsampled.
#' @export
#'
downSample <- function(qseaSet, nReads, renormalise = TRUE){
  counts <- qseaSet@count_matrix

  if (min(colSums(counts)) < nReads) {
    stop(glue::glue("Number of reads requested is less than the minimum {min(colSums(counts))}."))
  }

  message(glue::glue("Downsampling all samples to {nReads} each."))

  newCounts <- purrr::map_dfc(colnames(counts), function(colname){
        vec <- counts[,colname]

        sample(rep(1:length(vec), vec), replace = FALSE, size = nReads) %>%
          table() %>%
          tibble::enframe(name = "window") %>%
          dplyr::mutate(window = as.integer(window)) %>%
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

#' Calculate beta values for the windows covering probes from a methylation array.
#'
#' This function returns a wide data frame with beta values corresponding to overlapping probes from a methylation array.
#' This is suitable for input into tools in other packages that are designed for tables of array beta values.
#' Currently only "Infinium450k" is recognised, and only for GRCh38.
#'
#' @param qseaSet qseaSet object to calculate approximate beta values.
#' @param arrayDetails Either a recognised string or a GRanges object with an ID column.
#' @return A wide data frame with each row being a probe, and each sample being a column (plus the ID column)
#' @examples
#' convertToArrayBetaTable(exampleTumourNormal, arrayDetails = "Infinium450k")
#' convertToArrayBetaTable(exampleTumourNormal, arrayDetails = mesa::hg38_450kArrayGR)
#' @export
#'
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

#' This function calculates the fraction of reads in a qseaSet object which are present in a GRanges object
#' @param qseaSet A qseaSet
#' @param windowsToConsider Which windows to consider, as a GRanges object
#' @param numCountsNeeded Minimum number of reads required to consider that window
#' @export
#'
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

#' This function generates a data frame to use for the annotation of heatmaps, using Individual sample names not group
#' @param qseaSet A qseaSet object
#' @param ... Which columns of the sampleTable to use in the annotation data frame
#' @export
#'
getAnnotationDataFrameIndividual <- function(qseaSet, ...){

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

#' This function generates a data frame to use for the annotation of heatmaps, when using groups
#' @param qseaSet A qseaSet object
#' @param ... Which columns of the sampleTable to use in the annotation data frame
#' @export
#'
getAnnotationDataFrame <- function(qseaSet, ...){

  if (!("valid_fragments" %in% colnames(qsea::getSampleTable(qseaSet)))) {
    qseaSet <- qseaSet %>%
      addLibraryInformation()
  }

  qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(total_fragments = mean(total_fragments),
           relH = mean(relH)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(group,...) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("group") %>%
    return()
}

#' This function removes the normMethod suffix, if present, from the sample (or sample group) column names in a data frame of normalised values
#' @param dataTable  A data frame of normalised values for a set of windows (rows) and samples (columns), e.g. from [getDataTable()].
#' @param normMethod A character giving the normalisation method used (e.g. "beta" or "nrpm").
#' @return The dataTable with the normMethod suffix "_{normMethod}" removed from sample (or sample group) column names.
#' @export
#'
removeNormMethodSuffix <- function(dataTable, normMethod) {
  dplyr::rename_with(dataTable, ~ stringr::str_remove(.x, glue::glue("_{normMethod}(_means)?$")))
}

#' This function is a modified version of the [qsea::getPCA()] function
#' @param qseaSet A qseaSet object.
#' @param dataTable A data frame of normalised values for a set of windows (rows) and samples (columns), e.g. from [getDataTable()]. It must have seqnames, start and end columns. Can also be a [GenomicRanges::GRanges()] object with normalised values in the metadata columns.
#' @param regionsToOverlap Optional. Only windows in x overlapping `regionsToOverlap` will be considered. A [GenomicRanges::GRanges()] object or a data frame which can be coerced into a [GenomicRanges::GRanges()] object.
#' @param normMethod What normalisation method to use. Typically a character giving name of predefined normalisation method (e.g. "beta" or "nrpm"). See [qsea::normMethod()].
#' @param minEnrichment Minimum number of reads for beta values to not give NA. Passed to [getDataTable()].
#' @param useGroupMeans Whether to average samples over the group column (i.e. combine replicates)
#' @param minDensity A minimum CpG density level to filter out windows with values lower than.
#' @param topVarNum Number of most variable windows to keep. Defaults to 1000. If value is NA, NULL, Inf or is at least as big as the available number of windows, then all available windows are used for the PCA. Can also be a vector, in which case PCA is performed for each component.
#' @param topVarSamples Samples to use to determine variability. Either NULL or NA (use all samples; default), a character vector of sample names or a regular expression to match sample names on. Can also be a list, in which case PCA is performed for each list component. If `length(topVarNum) > 1`, list should be same length as `topVarNum` and each component will be matched with corresponding component of `topVarNum`.
#' @param center A logical value indicating if windows should be centred to have mean zero. Default is TRUE.
#' @param scale A logical value indicating if windows should be scaled to have unit variance. Default is FALSE.
#' @param nPC Number of principal components to be calculated. Default is 5.
#' @param returnDataTable A logical value indicating if the table of normalised values (generated by [getDataTable()] or by processing the input `dataTable`) should be returned as part of the output. See below for details.
#' @return A list with the following components:
#' \item{pcas}{A named list of objects of class prcomp. List names correspond to entries in the `params$topVar` component of the output.}
#' \item{samples}{Names of the samples analysed.}
#' \item{windows}{A list of names of the windows analysed. List names correspond to entries in the `params$topVar` component of the output.}
#' \item{dataTable}{If `returnDataTable = TRUE`, the table of normalised values. The table will contain windows after filtering based on `regionsToOverlap` and `minDensity`, but prior to filtering based on `topVarNum`. Windows with missing values will also have been removed. The table will also include column(s) of window standard deviations, where applicable, with column name(s) corresponding to entries in the `params$topVar` component of the output.}
#' \item{params}{A list of the parameters used.}
#' \item{windowFiltering}{A list of the initial number of windows and the number filtered out at each filtering step (not including filtering based on `topVarNum`).}
#' @export
#'
getPCA <- function(qseaSet, 
                   dataTable = NULL, 
                   regionsToOverlap = NULL, 
                   normMethod = "beta", 
                   minEnrichment = 3,
                   useGroupMeans = FALSE,
                   minDensity = 0,
                   topVarNum = 1000,
                   topVarSamples = NULL,
                   center = TRUE,
                   scale = FALSE,
                   nPC = 5,
                   returnDataTable = FALSE) {
  
  if (useGroupMeans) {
    groupString <- "sample group"
  } else {
    groupString <- "sample"
  }
  
  if (!is.null(dataTable)) {
    dataTable <- asValidGranges(dataTable)
    
    samples <- dataTable %>% 
      tidyr::as_tibble() %>% 
      dplyr::select(-tidyselect::any_of(c("seqnames", "start", "end", "width", "strand", "CpG_density"))) %>% 
      colnames()
    
    normMethodSuffixDetected <- stringr::str_detect(samples, glue::glue("_{normMethod}$"))
    
    if (all(normMethodSuffixDetected)) {
      samples <- stringr::str_remove(samples, glue::glue("_{normMethod}$"))
      dataTable <- removeNormMethodSuffix(dataTable %>% tidyr::as_tibble(), 
                                          normMethod) %>% 
        asValidGranges()
    } else if (any(normMethodSuffixDetected)) {
      stop(glue::glue("normMethod suffix '_{normMethod}' is present for some but not all of the {groupString} column names in dataTable."))
    }
    
    if (useGroupMeans) {
      if (!all(samples %in% names(qsea::getSampleGroups(qseaSet)))) {
        stop(glue::glue("At least one {groupString} name in dataTable does not have a matching {groupString} name in qseaSet.
             (If there is a normMethod suffix on the {groupString} column names in dataTable, check it matches the input normMethod argument: '_{normMethod}')"))
      }
    } else {
      if (!all(samples %in% qsea::getSampleNames(qseaSet))) {
        stop(glue::glue("At least one {groupString} name in dataTable does not have a matching {groupString} name in qseaSet.
             (If there is a normMethod suffix on the {groupString} column names in dataTable, check it matches the input normMethod argument: '_{normMethod}')"))
      }
    }
    
    if (length(setdiff(dataTable, getWindows(qseaSet))) > 0) {
      stop("At least one window in dataTable does not have a matching window in qseaSet.")
    }
    
    testSamples <- samples[1:min(5, length(samples))]
    inputValues <- dataTable[1:min(50, length(dataTable)), ]
    testValues <- getDataTable(qseaSet %>% 
                                 filter(sample_name %in% testSamples) %>% 
                                 filterByOverlaps(inputValues),
                               normMethod, 
                               useGroupMeans = useGroupMeans) %>% 
      dplyr::arrange(seqnames, start, end)
    
    if (!isTRUE(all.equal(testValues, inputValues %>% 
                          tidyr::as_tibble() %>% 
                          dplyr::select(tidyselect::all_of(colnames(testValues))) %>% 
                          dplyr::arrange(seqnames, start, end)))) {
      warning("Newly-generated dataTable from qseaSet on a small subset of samples/windows does not match input dataTable.
              \n")
    }
    initialNumWindows <- length(dataTable)
  } else {
    
    if (useGroupMeans) {
      samples <- names(qsea::getSampleGroups(qseaSet))
    } else {
      samples <- qsea::getSampleNames(qseaSet)
    }
    initialNumWindows <- getWindows(qseaSet) %>% length()
  }
  
  message(glue::glue("------------------------------
                     Initial number of windows = {initialNumWindows}."))
  
  if (!is.null(regionsToOverlap)) {
    regionsToOverlap <- regionsToOverlap %>% 
      tibble::as_tibble() %>% 
      dplyr::select(tidyselect::any_of(c("seqnames", "start", "end", "CpG_density"))) %>%  # keep only minimum columns necessary
      asValidGranges()
    
    if (is.null(dataTable)) {
      qseaSet <- qseaSet %>% 
        filterByOverlaps(regionsToOverlap = regionsToOverlap)
      
      numWindowsRemovedRegionOverlap <- initialNumWindows - length(getWindows(qseaSet))
      message(glue::glue("Filtered out {numWindowsRemovedRegionOverlap} windows using regionsToOverlap: {length(getRegions(qseaSet))} windows remaining."))
      
    } else {
      dataTable <- dataTable %>% 
        plyranges::filter_by_overlaps(y = regionsToOverlap)
      
      numWindowsRemovedRegionOverlap <- initialNumWindows - length(dataTable)
      message(glue::glue("Filtered out {numWindowsRemovedRegionOverlap} windows using regionsToOverlap: {length(dataTable)} windows remaining."))
    }
  } else {
    numWindowsRemovedRegionOverlap <- NULL
  }
  
  if (is.null(dataTable)) {
    message("-----------")
    dataTable <- qseaSet %>%
      filterWindows(CpG_density >= minDensity) %>%
      getDataTable(normMethod = normMethod, useGroupMeans = useGroupMeans)
    message("-----------")
    
    if (minDensity > 0) {
      numWindowsRemovedMinDensity <- length(getWindows(qseaSet)) - nrow(dataTable)
      message(glue::glue("Filtered out {numWindowsRemovedMinDensity} windows with CpG_density < {minDensity}: {nrow(dataTable)} windows remaining."))
    } else {
      numWindowsRemovedMinDensity <- NULL
    }
    
  } else {
    currentNumWindows <- length(dataTable)
    dataTable <- dataTable %>% 
      tidyr::as_tibble() %>% 
      dplyr::left_join(qseaSet %>% 
                         getWindows() %>% 
                         tidyr::as_tibble() %>% 
                         dplyr::select(seqnames, start, end, CpG_density)) %>% 
      dplyr::filter(CpG_density >= minDensity)
    
    if (minDensity > 0) {
      numWindowsRemovedMinDensity <- currentNumWindows - nrow(dataTable)
      message(glue::glue("Filtered out {numWindowsRemovedMinDensity} windows with CpG_density < {minDensity}: {nrow(dataTable)} windows remaining."))
    } else {
      numWindowsRemovedMinDensity <- NULL
    }
    
  }
  
  currentNumWindows <- nrow(dataTable)
  
  dataTable <- dataTable %>% 
    tidyr::drop_na(tidyr::all_of(samples))
  
  numWindowsRemovedMissingVals <- currentNumWindows - nrow(dataTable)
  
  if (numWindowsRemovedMissingVals > 0) {
    message(glue::glue("Filtered out {numWindowsRemovedMissingVals} windows with at least one missing value: {nrow(dataTable)} windows remaining.\n
                       ------------------------------
                       \n"))
  } else {
    message(glue::glue("No windows have missing values.\n
            ------------------------------
            \n"))
  }  

  
  if (is.null(topVarNum)) {
    topVarNum <- NA
  } else if (!(length(topVarNum) == 1 && is.na(topVarNum)) & !is.vector(topVarNum, mode = "numeric")) {
    stop("topVarNum should be a numeric vector (or NULL")    
  }
  
  if (!is.list(topVarSamples)) {
    topVarSamples <- list(topVarSamples)
  }
  
  if (length(topVarNum) > 1 && !(length(topVarSamples) %in% c(1, length(topVarNum)))) {
    stop("If topVarSamples is a list and length(topVarNum) > 1, topVarSamples should be the same length as topVarNum.")
  }
  
  
  # replace NA with NULL
  topVarSamples <- purrr::map(topVarSamples, function(tVS) {
    
    if (length(tVS) == 1 && is.na(tVS)) {
      tVS <- NULL
    }
    
    return(tVS)
  })
  
  topVarSamplesInput <- topVarSamples
  
  topVarSamples <- purrr::map(topVarSamples, function(tVS) {
    
    if (!is.null(tVS)) {
      if (!is.vector(tVS, mode = "character")) {
        stop(glue::glue("topVarSamples should be NULL, a character vector of {groupString} names or a regular expression (or a list of these)."))
        
      } else if (length(tVS) > 1) { # character vector of sample (group) names
        notInSamples <- setdiff(tVS, samples)
        if (length(notInSamples) > 0) {
          stop(glue::glue("topVarSamples contains {groupString} names that are not in the qseaSet and/or dataTable:
                        {paste0(notInSamples, collapse = ', ')}."))
        }
        
      } else { # regular expression to match
        tVS <- stringr::str_subset(samples, tVS)
        
      }
      
    } else {
      tVS <- samples
      
    }
    
    return(tVS)
    
  })
  
  topVar <- tibble::tibble(topVarNum = topVarNum, topVarSamples = topVarSamples, topVarNumInput = topVarNum, topVarSamplesInput = topVarSamplesInput) %>% 
    dplyr::mutate(topVarNum = ifelse(topVarNum >= nrow(dataTable), NA, topVarNum),
                  topVarSamples = ifelse(is.na(topVarNum), list(NULL), topVarSamples),
                  inputChanged = !purrr::map2_lgl(topVarNum, topVarNumInput, identical) | 
                    !purrr::map2_lgl(topVarSamples, topVarSamplesInput, identical))
  
  if (any(topVar$topVarNumInput > nrow(dataTable), na.rm = TRUE)) {
    message(glue::glue("The following topVarNum values are larger than the number of remaining windows (= {nrow(dataTable)}): {paste0(dplyr::filter(topVar, topVarNumInput > nrow(dataTable)) %>% pull(topVarNumInput) %>% unique(), collapse = ', ')}"))
    
    if (any(is.na(topVar$topVarNumInput) | topVar$topVarNumInput == nrow(dataTable))) {
      message("These values are not used; PCA is already being done with all remaining windows.\n")
    } else {
      message("These values are not used; PCA will be done with all remaining windows instead.\n")
    }
  }
  
  
  topVar <- topVar %>% 
    dplyr::distinct(topVarNum, topVarSamples, .keep_all = TRUE) %>% 
    dplyr::group_by(topVarSamples) %>% 
    dplyr::mutate(windowSdName = ifelse(!is.na(topVarNum), glue::glue("windowSd{dplyr::cur_group_id()}"), NA), .before = 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(windowSdName, topVarNum) %>% 
    dplyr::mutate(pcaName = glue::glue("pca{dplyr::row_number()}"), .before = 1) %>% 
    dplyr::bind_rows(topVar) %>% 
    dplyr::distinct(topVarNum, topVarSamples, topVarNumInput, topVarSamplesInput, .keep_all = TRUE)
  
  rowSds <- function(x) {
    sqrt(rowSums((x - rowMeans(x)) ^ 2) / (ncol(x) - 1))
  }
  
  if (any(!is.na(topVar$windowSdName))) {
    dataTable <- topVar %>% 
      tidyr::drop_na(windowSdName) %>% 
      dplyr::select(topVarSamples, windowSdName) %>% 
      dplyr::distinct() %>%
      dplyr::mutate(topVarSamples = purrr::set_names(topVarSamples, windowSdName)) %>% 
      dplyr::pull(topVarSamples) %>% 
      purrr::imap(function(tVS, nm) {
        
        if (length(setdiff(samples, tVS)) == 0) {
          message(glue::glue("Calculating standard deviation for each window across all {length(samples)} {groupString}s:
                           {paste0(tVS, collapse = ', ')}.
                           -> column name {nm}.
                           \n"))
        } else {
          message(glue::glue("Calculating standard deviation for each window across {length(tVS)} of {length(samples)} {groupString}s:
                           {paste0(tVS, collapse = ', ')}.
                           -> column name {nm}.
                           \n"))
        }
        
        if (length(tVS) <= 2) {
          stop(glue::glue("Standard deviation calculated on less than 3 {groupString} names. Insufficent number of {groupString}s to calculate variance."))
        }
        
        dataTable <- dataTable %>%  
          dplyr::mutate({{nm}} := dplyr::select(., tidyr::all_of(tVS)) %>% 
                          rowSds()) %>% 
          dplyr::select(seqnames, start, end, tidyr::all_of(nm))
        
      }) %>%
      purrr::reduce(dplyr::left_join) %>% 
      dplyr::left_join(dataTable, .)
  }
  
  pca <- topVar %>% 
    dplyr::filter(!is.na(pcaName)) %>% 
    dplyr::group_by(windowSdName) %>% 
    dplyr::group_map(.keep = TRUE, .f = function(sdGp, sdName) {
      
      if (!is.na(sdName$windowSdName)) {
        dataTable <- dataTable %>%  
          dplyr::arrange(dplyr::desc(!!dplyr::sym(sdName$windowSdName)))
      }
      
      purrr::pmap(sdGp, function(pcaName, windowSdName, topVarNum, topVarSamples, ...) {
        
        if (!is.na(windowSdName)) {
          th <- dataTable %>% 
            dplyr::pull({{windowSdName}}) %>% 
            dplyr::nth(topVarNum)
          
          dataTable <- dataTable %>% 
            dplyr::filter(!!dplyr::sym(windowSdName) >= th)
          
          if (length(topVarSamples) == length(samples)) {
            
            message(glue::glue("------------------------------\n
                             Filtering windows based on standard deviation across all {length(topVarSamples)} {groupString}s ({windowSdName}). 
                             Standard deviation threshold = {format(th, digits = 3)} resulting in {nrow(dataTable)} windows.
                             \n"))
          } else {
            message(glue::glue("------------------------------\n
                             Filtering windows based on standard deviation across {length(topVarSamples)} {groupString}s ({windowSdName}). 
                             Standard deviation threshold = {format(th, digits = 3)} resulting in {nrow(dataTable)} windows.
                             \n"))
          }
          
          
        } else {
          
          message(glue::glue("------------------------------\n
                             No filtering of windows based on window standard deviation.
                             \n"))
          th <- NA
        }
        
        dataTable <- dataTable %>% 
          dplyr::mutate(window = getWindowNames(.)) %>% 
          tibble::column_to_rownames("window") %>% 
          dplyr::select(tidyr::all_of(samples))
        
        message(glue::glue("Performing PCA with {ncol(dataTable)} {groupString}s and {nrow(dataTable)} windows
                           -> {pcaName}.
                             \n"))
        
        dataTable <- dataTable %>% 
          t()
        
        if (!is.numeric(dataTable)) {
          stop("Input to PCA contains non-numeric values.")
        }
        
        if (any(is.na(dataTable))) {
          stop("Input to PCA contains missing values.")
        }
        
        if (any(!is.finite(dataTable))) {
          stop("Input to PCA contains infinite values.s")
        }
        
        prcompObj <- dataTable %>% 
          stats::prcomp(center = center, scale. = scale, rank. = nPC)
        
        return(list(prcompObj = prcompObj, th = th))
        
      }) %>% 
        purrr::set_names(sdGp$pcaName)
    }) %>% 
    purrr::list_flatten()
        
  message("------------------------------")
  
  th <- purrr::map(pca, "th") %>% 
    unlist()
  
  pcas <- purrr::map(pca, "prcompObj")
  
  windows <- purrr::map(pcas, ~ rownames(.x$rotation))
    
  return(list(pcas = pcas, 
              samples = samples,
              windows = windows,
              dataTable = if (returnDataTable) {
                dataTable
                } else {
                  NULL
                }, 
              params = list(regionsToOverlap = regionsToOverlap,
                            normMethod = normMethod,
                            minEnrichment = minEnrichment,
                            useGroupMeans = useGroupMeans,
                            minDensity = minDensity,
                            topVar = topVar,
                            windowSdThreshold = th,
                            center = center,
                            scale = scale,
                            nPC = nPC),
              windowFiltering = list(initial = initialNumWindows,
                                     notInRegionsToOverlap = numWindowsRemovedRegionOverlap,
                                     belowMinDensity = numWindowsRemovedMinDensity,
                                     containMissingVals = numWindowsRemovedMissingVals))) 
}


#' This function takes a qseaSet and sums over nrpm within a GRanges object
#' @param qseaSet The qseaSet object.
#' @param GRanges A GRanges object to sum over
#' @param samples What samples to use
#' @param normMethod What normalisation method to use
#' @param cutoff What cutoff to use to count the windows above
#' @return A table
#' @export
#'
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


#' This function takes a qseaSet and makes a wide transposed table, suitable for using for machine learning
#' @param qseaSet The qseaSet object.
#' @param normMethod What normalisation method to use
#' @param ... Other columns to add from the sampleTable in the output (along with sample_name). Must specify normMethod if using. Can use everything() for all columns.
#' @return A table
#' @export
#'
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


#' This function takes a qseaSet and makes a table of nrpm values
#' @param qseaSet The qseaSet object.
#' @param useGroupMeans Whether to give means of the group column, rather than individual samples.
#' @return A table of data, one row per window
#' @export
#'
getCountTable <- function(qseaSet, useGroupMeans = FALSE){

  if(useGroupMeans){

    qseaSet %>%
      qsea::makeTable(groupMeans =  qsea::getSampleGroups(.), norm_methods = "counts") %>%
      dplyr::rename_with(~ stringr::str_replace_all(.x, "_counts_means", "")) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble() %>%
      return()

  } else {

    qseaSet %>%
      qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = "counts") %>%
      dplyr::rename_with(~ stringr::str_replace_all(.x, "_counts", "")) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble() %>%
      return()
  }

}


#' This function takes a qseaSet and makes a table of nrpm values
#' @param qseaSet The qseaSet object.
#' @param useGroupMeans Whether to give means of the group column, rather than individual samples.
#' @return A table of data, one row per window
#' @export
#'
getNRPMTable <- function(qseaSet, useGroupMeans = FALSE){

  if(useGroupMeans){

    qseaSet %>%
      qsea::makeTable(groupMeans =  qsea::getSampleGroups(.), norm_methods = "nrpm") %>%
      dplyr::rename_with(~ stringr::str_replace_all(.x, "_nrpm_means", "")) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble() %>%
      return()

  } else {

    qseaSet %>%
      qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = "nrpm") %>%
      dplyr::rename_with(~ stringr::str_replace_all(.x, "_nrpm", "")) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble() %>%
      return()
  }

}

#' This function takes a qseaSet and makes a table of beta values
#' @param qseaSet The qseaSet object.
#' @param useGroupMeans Whether to give means of the group column, rather than individual samples.
#' @return A table of data, one row per window
#' @export
#'
getBetaTable <- function(qseaSet, useGroupMeans = FALSE){

  if(useGroupMeans){

    qseaSet %>%
      qsea::makeTable(groupMeans =  qsea::getSampleGroups(.), norm_methods = "beta") %>%
      dplyr::rename_with(~ stringr::str_replace_all(.x, "_beta_means", "")) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble() %>%
      return()

  } else {

    qseaSet %>%
      qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = "beta") %>%
      dplyr::rename_with(~ stringr::str_replace_all(.x, "_beta", "")) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble() %>%
      return()
  }

}

#' This function takes a qseaSet and calculates a summary statistic across the windows.
#' @param qseaSet The qseaSet object.
#' @param regionsToOverlap A GRanges object (or dataframe coercible to one) containing the windows to summarise over. If missing, will use all the regions in the qseaSet.
#' @param fn A function to apply across the windows. e.g. mean, median, sd.
#' @param suffix A suffix for adding to the new columns of data, to clarify what regions for instance.
#' @param addSampleTable A boolean with whether to add the sampleTable on to the output table
#' @param normMethod One or more normalisation methods to use, e.g. "nrpm" or "beta" or c("nrpm", "beta").
#' @param naMethod What method to use to deal with NA values (for beta values). Options are "drop" or "na.rm" (default).
#' @param minEnrichment For beta values only, the minimum number of reads required for a window to be fully methylated, qsea replaces with NA below this value.
#' @param fnName Name of the function. Should only be necessary to use if you are doing something unusual, otherwise detected automatically.
#' @return A table of data, one row per sample, with columns indicating the application of the summary statistic to the
#' @export
#'
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
      qsea::makeTable(norm_methods = normMethod,
                      samples = qsea::getSampleNames(qseaSet),
                      minEnrichment = minEnrichment)

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

#' This function takes a qseaSet and adds to the sampleTable summary statistics calculated over a set of windows.
#' @param qseaSet The qseaSet object.
#' @param regionsToOverlap A GRanges object (or data frame coercible to one) containing the windows to summarise over. If missing, will use all the regions in the qseaSet.
#' @param fn A function to apply across the windows. e.g. mean, median, sd.
#' @param suffix A suffix for adding to the new columns of data, to clarify where the regions came from for instance.
#' @param normMethod One or more normalisation methods to use, e.g. "nrpm" or "beta" or c("nrpm", "beta").
#' @param naMethod What method to use to deal with NA values (for beta values). Options are "drop" or "impute".
#' @param minEnrichment For beta values only, the minimum number of reads required for a window to be fully methylated, qsea replaces with NA below this value.
#' @return A table of data, one row per sample, with columns indicating the application of the summary statistic to the
#' @export
#'
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
    left_join(summaryTable, by = sample_name)

  return(qseaSet)

}
#' This function takes a qseaSet and calculates some stats as to which genomic regions the reads lie.
#' @param qseaSet The qseaSet object.
#' @param cutoff The value required to call a window as being above that cutoff or not
#' @param normMethod What normalisation method to use (nrpm or beta)
#' @param minEnrichment Number of reads required for a fully methylated region to not be masked (put NA) beta values only.
#' @return A table of summarised data, with windows split out into both the landscape (island/shore/shelf) and genomic region (promoter/exon/intron etc)
#' @export
#'
getGenomicFeatureDistribution <- function(qseaSet, cutoff = 1 , normMethod = "nrpm", minEnrichment = 3){
  temp <- qseaSet %>%
    qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = normMethod, minEnrichment = minEnrichment) %>%
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

#' This function takes a qseaSet and generates a table of data containing data for each sample
#' @param qseaSet The qseaSet object.
#' @param normMethod What normalisation method to use (nrpm or beta), can be given multiple.
#' @param useGroupMeans Whether to use the group column to average over replicates.
#' @param minEnrichment Minimum number of reads for beta values to not give NA
#' @param addMethodSuffix Whether to include a suffix corresponding to the normalisation method, such as Sample1_beta. This suffix is always present if multiple normalisationMethods are given.
#' @return A table of data, with a column for each individual sample
#' @export
#'
getDataTable <- function(qseaSet, normMethod = "nrpm", useGroupMeans = FALSE, minEnrichment = 3, addMethodSuffix = FALSE){

  if(useGroupMeans){

    tab <- qseaSet %>%
      qsea::makeTable(groupMeans =  qsea::getSampleGroups(.), norm_methods = normMethod, minEnrichment = minEnrichment) %>%
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

    tab <- qseaSet %>%
      qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = normMethod, minEnrichment = minEnrichment) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble()

    if(length(normMethod) == 1 & !addMethodSuffix) {
      tab <- tab %>%
        dplyr::rename_with(~ stringr::str_replace_all(.x, glue::glue("_{normMethod}"), ""))
    }

    return(tab)
  }

}

#' This function takes a qseaSet and writes individual bigWig files with the coverage over each window.
#' @param qseaSet The qseaSet object.
#' @param folderName A folder to save the output bigWig files into
#' @param normMethod What normalisation method to use (e.g. nrpm or beta)
#' @param useGroupMeans Whether to average over the replicates in the group
#' @param naVal A value to replace NA values with (for beta values).
#' @return A set of bigWig files for each sample in the qseaSet, with coverage over all the windows in the qseaSet
#' @export
#'
writeBigWigs <- function(qseaSet, folderName, normMethod = "nrpm", useGroupMeans = FALSE, naVal = -1){

  dir.create(folderName, showWarnings = FALSE)

  dataTable <- qseaSet %>%
    getDataTable(normMethod = normMethod, useGroupMeans = useGroupMeans)

  if(!useGroupMeans) {
    mapNames <- qseaSet %>%
      qsea::getSampleNames()

  } else {
    mapNames <- qseaSet %>%
      qsea::getSampleGroups()

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

#' This function takes a qseaSet and sets all the library factor columns to one, including in the sampleTable.
#' @param qseaSet The qseaSet object.
#' @return A table of data, with a column for each individual sample
#' @export
#'
removeLibraryFactors <- function(qseaSet){
  qseaSet <- qseaSet %>%
    qsea::addLibraryFactors(1)

  if("library_factor" %in% colnames(qsea::getSampleTable(qseaSet))){
    qseaSet <- qseaSet %>%
      dplyr::mutate(library_factor = 1)
  }
return(qseaSet)
}

