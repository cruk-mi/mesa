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
#' This function extracts the information from one contrast into a wide table
#'
#' @param dataTable The output of makeTable
#' @param genome A genome string to set the rest of the parameters (currently only hg38/GRCh38 supported)
#' @param TxDb A TxDb database object (unquoted) to pass to ChIPseeker::annotatePeak
#' @param annoDb A string giving a Bioconductor annotation package, such as "org.Hs.eg.db"
#' @param CpGislandsGR A GRanges object giving locations of CpG islands
#' @param FantomRegionsGR A GRanges object giving Fantom enhancer regions
#' @return A tibble with the data, augmented with ChIPseeker region location and CpG island information.
#' @export
annotateData <- function(dataTable, genome = "hg38", TxDb = NULL, annoDb = NULL, CpGislandsGR = NULL,
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

    annoDb = "org.Hs.eg.db" }

  if(genome  %in% c("hg38","GRCh38") & is.null(CpGislandsGR)) { CpGislandsGR = mesa::hg38CpGIslands }

  if(genome  %in% c("hg38","GRCh38") & is.null(FantomRegionsGR)) { FantomRegionsGR = mesa::FantomRegions %>% plyranges::as_granges()}

  chipseekerData <- dataTable %>%
    qseaTableToChrGRanges() %>%
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

#' This function takes a qseaSet and filters out windows which have expression in a defined set of samples
#'
#' @param qseaSet The qseaSet object.
#' @param samplesToFilterOut A set of sample names to filter out, or a string to match in the sample names
#' @param maxValue The maximum value to allow in any of the samples
#' @param normMethod The type of normalisation method to use (nrpm, beta, count)
#' @param .swap Whether to change the printed messages as we are actually going to be keeping these windows rather than filtering them.
#' @return A qseaSet object, with only a subset of the windows.
#' @export
removeWindowsOverCutoff <- function(qseaSet, samplesToFilterOut, maxValue, normMethod, .swap = FALSE){

  rowAny <- function(x) {rowSums(x) > 0}

  if (!(normMethod %in% c("beta","nrpm","counts"))) {
    stop(glue::glue("normMethod must be either beta, nrpm or counts, not {normMethod}"))
  }

  ##TODO Make work for beta values
  if (normMethod == "beta") {
    stop(glue::glue("normMethod of beta doesn't currently work, need to deal with NAs."))
  }

  samplesNotInSet <- setdiff(samplesToFilterOut, qsea::getSampleNames(qseaSet))

  if (length(samplesToFilterOut) == 1 & length(samplesNotInSet) > 0 & is.character(samplesToFilterOut)) {
    sampleNameString <- samplesToFilterOut
    samplesToFilterOut <- stringr::str_subset(qsea::getSampleNames(qseaSet), samplesToFilterOut)
    message(glue::glue("Considering {length(samplesToFilterOut)} samples containing {sampleNameString} in the name."))
  } else if (length(samplesNotInSet) > 0 ) {
    stop(glue::glue("Sample {samplesNotInSet} not present in the qseaSet!

                    "))
  }

  groupNames <- qsea::getSampleGroups(qseaSet)[unlist(lapply(qsea::getSampleGroups(qseaSet), function(x) any(samplesToFilterOut %in% x)))]

  dataTable <- suppressMessages(qsea::makeTable(qseaSet, norm_methods = normMethod,
                               groupMeans = groupNames
    ))

  keep <- dataTable %>%
    as.data.frame() %>%
    tibble::rowid_to_column(var = ".rowID") %>%
    dplyr::select(tidyselect::matches("_means")) %>%
    apply(1,max) %>%
    {which(. <= maxValue)}

  if(!.swap){
     message(glue::glue("Removing {nrow(dataTable) - length(keep)} windows based on {length(samplesToFilterOut)} samples, {length(keep)} remaining"))
  }else{
    message(glue::glue("Keeping {nrow(dataTable) - length(keep)} windows based on {length(samplesToFilterOut)} samples, {length(keep)} removed"))
  }
  return(invisible(filterByOverlaps(qseaSet, keep)))

}


#' This function takes a qseaSet and keeps the windows which have expression above a cutoff in a defined set of samples
#'
#' @param qseaSet The qseaSet object.
#' @param samplesToFilterOn A set of sample names to keep, or a string to match in the sample names
#' @param minValue The minimum value to require in any of the samples
#' @param normMethod The type of normalisation method to use (nrpm, beta, count)
#' @return A qseaSet object, with only a subset of the windows.
#' @export
keepWindowsOverCutoff <- function(qseaSet, samplesToFilterOn, minValue, normMethod){

  # find which regions would be removed by removeWindowsOverCutoff and keep the opposite set!
  qseaSet %>%
    removeWindowsOverCutoff(samplesToFilterOut = samplesToFilterOn, maxValue = minValue, normMethod = normMethod, .swap = TRUE) %>%
    qsea::getRegions() %>%
    filterByNonOverlaps(qseaSet, .) %>%
    return()
}


#' This function takes a qseaSet and finds which windows have more reads than would be expected if the reads followed a Poisson distribution over the whole genome.
#'
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

#' This function takes a qseaSet and makes a new sample by mixing two samples
#' @param qseaSet The qseaSet object.
#' @param arrayReadTable Data frame with the array probe beta values
#' @param arraySample Name of the array sample to mix
#' @param qseaSample Name of the qseaSet sample to mix with
#' @param nReadsTotal Number of reads in total to have after mixing
#' @param proportion The proportion to take from sample1, the rest will come from sample2
#' @param newName A name to give the new sample
#' @param groupName A name to use in the group column in the sampleTable
#' @param onlyNew Whether to only return the new sample.
#' @return A qseaSet object with an extra
#' @export
#'

mixArrayWithQset <- function(qseaSet, arrayReadTable, arraySample, qseaSample, nReadsTotal, proportion,
                             newName = NULL,
                             groupName = NULL,
                             onlyNew = FALSE){

  if(!(arraySample %in% colnames(arrayReadTable))){
    stop(glue::glue("Column {arraySample} not present in the array Table"))
  }

  qsea:::checkSamples(qseaSet,c(qseaSample))

  if(is.null(newName)) { newName <- paste0("Mix","_",arraySample,"_",qseaSample,"_",proportion) }

  qseaSetFiltered <- qseaSet %>%
    filterByOverlaps(arrayReadTable %>% tidyr::drop_na())

  message(newName)

  if(is.null(groupName)) { groupName <- newName }

  nReads1 <- ceiling(nReadsTotal * proportion)
  nReads2 <- nReadsTotal - nReads1

  if (proportion < 0 || proportion > 1) {
    stop(glue::glue("Proportion must be a number between 0 and 1, not {proportion}."))
  }

  validFragNums <- qseaSetFiltered@libraries$file_name[c( qseaSample),"valid_fragments", drop = FALSE] %>%
    tibble::rownames_to_column() %>%
    tibble::deframe()

  counts1 <- arrayReadTable %>%
    tidyr::drop_na() %>%
    dplyr::select(tidyselect::matches(paste0("^",arraySample,"$"))) %>%
    as.data.frame()

  colnames(counts1) <- arraySample

  counts2 <- qseaSetFiltered@count_matrix[,c(qseaSample), drop = FALSE]

  counts <- cbind(counts1,counts2)

  sampleColSums <- colSums(counts)

  message(glue::glue("Mixing {nReads1} reads from {arraySample} with {nReads2} reads from {qseaSample}."))

  onTargetFrac <- sampleColSums[qseaSample]/validFragNums[qseaSample]

  if (sampleColSums[arraySample]  < nReads1 * onTargetFrac) {
    stop(glue::glue("Not enough reads in {arraySample}, only {sampleColSums[arraySample]} out of {nReads1 * onTargetFrac} requested."))
  }

  if (sampleColSums[qseaSample]  < nReads2 * onTargetFrac) {
    stop(glue::glue("Not enough reads in {qseaSample}, only {sampleColSums[qseaSample]} out of {nReads2 * onTargetFrac} requested."))
  }

  windows1 <- sample(rep(1:nrow(counts), counts[,arraySample]), replace = FALSE, size = nReads1 * onTargetFrac)
  windows2 <- sample(rep(1:nrow(counts), counts[,qseaSample]), replace = FALSE, size = nReads2 * onTargetFrac)

  newCounts <- c(windows1, windows2) %>%
    tibble::enframe(name = NULL, value = "window") %>%
    dplyr::count(window, name = "n") %>%
    dplyr::left_join(tibble::tibble(window = 1:nrow(counts)),., copy = TRUE, by = "window") %>%
    dplyr::mutate(n = tidyr::replace_na(n,0)) %>%
    dplyr::select(n)

  newSet <- qseaSetFiltered %>%
    subsetQset(samplesToKeep = qsea::getSampleNames(qseaSetFiltered)[1]) %>%
    renameQsetNames(paste0("^",qsea::getSampleNames(qseaSetFiltered)[1],"$"), newName)

  newSet@count_matrix <- as.matrix(newCounts)
  colnames(newSet@count_matrix) <- newName

  reducedCNVMat <- qseaSet %>%
    qsea::getCNV() %>%
    GenomicRanges::mcols() %>%
    as.matrix()

  GenomicRanges::mcols(newSet@cnv) <- reducedCNVMat[,qseaSample]
  colnames(GenomicRanges::mcols(newSet@cnv)) <- newName

  newSet@zygosity <- qseaSet@zygosity[qseaSample,,drop = FALSE]
  rownames(newSet@zygosity) <- newName

  newSet@sampleTable <- tibble::tibble(sample_name = newName,
                                       rownameCol = newName,
                                       group = groupName,
                                       #type = qseaSet %>% qsea::getSampleTable() %>% dplyr::filter(sample_name == !!arraySample) %>% dplyr::pull(type),
                                       #tumour = qseaSet %>% qsea::getSampleTable() %>% dplyr::filter(sample_name == !!arraySample) %>% dplyr::pull(tumour),
                                       sample1 = arraySample,
                                       sample2 = qseaSample,
                                       prop1 = proportion) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("rownameCol") %>%
    as.data.frame()

  newSet@libraries$file_name[,"valid_fragments"] <- nReadsTotal
  newSet@libraries$file_name[,"offset"] <- NA
  newSet@libraries$file_name[,"library_factor"] <- NA

  rownames(newSet@libraries$file_name) <- newName

  # if(renormalise){
  #   newSet <- addQseaNormalisationSteps(newSet, enrichmentMethod = "blind1-15")
  # }

  if (onlyNew) {
    return(newSet)
  } else {
    return(combineQsets(qseaSetFiltered, newSet))
  }

}





#' This function takes a qseaSet and downsamples the reads.
#' @param qseaSet The qseaSet object.
#' @param nReads How many reads to downsample to.
#' @param renormalise Whether to renormalise the resulting qseaSet.
#' @return A qseaSet object with the reads downsampled.
#' @export
#'
downsampleQsea <- function(qseaSet, nReads, renormalise = TRUE){
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

  if (renormalise) {   qseaSet <- addQseaNormalisationSteps(qseaSet)  }

  return(qseaSet)
}

#' This function takes a qseaSet and sums over nrpm within a GRanges object
#' @param qseaSet The qseaSet object.
#' @param GRanges A GRanges object to sum over
#' @param samples What samples to use
#' @param subtractLevel A value to subtract from the nrpm for each window before summing over them. Maybe useful to remove background noise.
#' @return A table
#' @export
#'
getNormalisedReadSum <- function(qseaSet, GRanges, samples = NULL, subtractLevel = 0){

  if (is.null(samples)) {
    samples <- qsea::getSampleNames(qseaSet)
  }

  clip <- function(x, a, b) {
    a + (x - a > 0) * (x - a) - (x - b > 0) * (x - b)
  }

  qseaSet %>%
    filterByOverlaps(GRanges) %>%
    qsea::makeTable(norm_methods = "nrpm", samples = samples) %>%
    dplyr::select(dplyr::matches("nrpm")) %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "_beta|_nrpm|_means")) %>%
    {clip(. - subtractLevel, a = 0, b = 1000000)} %>%
    colSums() %>%
    tibble::enframe(name = "sample_name", value = "nrpmSum") %>%
    dplyr::left_join(qsea::getSampleTable(qseaSet))

}



#' This function calculates the beta values for the windows covering the hg38 array probes.
#' @param qseaSet qseaSet object to calculate the probe beta values for
#' @export
#'
calculateArrayBetas <- function(qseaSet) {

  qseaSet %>%
    qsea::makeTable(norm_methods = "beta", samples = qsea::getSampleNames(.), ROIs = mesa::hg38_450kArrayGR) %>%
    dplyr::select(-tidyselect::matches("ROI_start|ROI_end|ROI_chr")) %>%
    dplyr::rename(ID = ROI_ID) %>%
    dplyr::select(ID, tidyselect::matches("beta")) %>%
    dplyr::select_if(function(x){!all(is.na(x))}) %>%
    dplyr::rename_with( ~ stringr::str_remove(., "_means"), tidyselect::matches("_means")) %>%
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



#' This function edits the getPCA function from qsea to allow for batch effect correction
#' @param qs A qseaSet object
#' @param chr Which chromosomes to use
#' @param minRowSum Minimum number of reads in a window across all the samples to keep it
#' @param minEnrichment Minimal number of expected reads for a fully methylated window (for transformation to absolute methylation level)
#' @param keep a vector of row indices to keep (optional)
#' @param normMethod What normalisation method to use
#' @param topVar Number of most variable windows to keep
#' @param samples Which samples to use
#' @param batchVariable A variable in the sampleTable to batch correct for, using limma::removeBatchEffect
#' @export
#'
getPCAwithBatch <- function(qs, chr = qsea::getChrNames(qs), minRowSum = 20, minEnrichment = 3, keep,
                            normMethod = "beta",
                            topVar = 1000, samples = qsea::getSampleNames(qs),
                            batchVariable = NULL){

  if( is.null(samples)) { samples =  qsea::getSampleNames(qs)}

  samples = qsea:::checkSamples(qs, samples)

  if( missing(keep)) {
    keep = which(rowSums(qsea::getCounts(qs)) >= minRowSum )
  }
  else {
    keep = intersect(keep, which(rowSums(qsea::getCounts(qs)) >= minRowSum ))
  }

  if (methods::is(normMethod, "character")) {
    normMethod = qsea::normMethod(normMethod)
  }

  if( !all( qsea::getChrNames(qs) %in% chr )){
    keep = keep[as.character(GenomeInfoDb::seqnames(qsea::getRegions(qs)[keep])) %in% chr]}

  if( length(keep) <= 10) {
    stop("not enough regions left:",length(keep)," Regions, minimum is 10")
  }

  vals = qsea:::getNormalizedValues(qs, methods = normMethod, minEnrichment = minEnrichment,
                                    windows = keep, samples = samples)

  missing = rowSums(is.na(vals)) > 0
  vals = vals[!missing,] - rowMeans(vals[!missing,])
  if ( any(!is.finite(vals))) {
    stop("Infinite values due to log or logit transformation. ",
         "Please specify minVal and maxVal.")
  }
  if ( !is.null(topVar) && nrow(vals) > topVar) {
    #rv=apply(vals,1,var) #this is slow
    #rv=rowSums((vals - rowMeans(vals))^2)/(dim(vals)[2] - 1)
    #since only the order matters:
    rv = rowSums((vals - rowMeans(vals))^2)
    th = sort(rv,decreasing = TRUE)[topVar]
    vals = vals[rv >= th,]
    keep = keep[rv >= th]
  }
  #make names for the selected regions
  names = paste0(GenomeInfoDb::seqnames(qsea::getRegions(qs)[keep]), ":",
                 IRanges::start(qsea::getRegions(qs)[keep]), "-", IRanges::end(qsea::getRegions(qs)[keep]))

  if ( !is.null(batchVariable)) {

    numBatches <- qs %>% pullQset(batchVariable) %>% unique()

    if (length(numBatches) > 1) {
      message(glue::glue("Normalising PCA based on {length(numBatches)} values of {batchVariable}"))
      vals <- vals %>% limma::removeBatchEffect(batch = dplyr::pull(qsea::getSampleTable(qs), batchVariable))
    }
  }

  svdVals = svd(vals)
  methods::new('qseaPCA', svd = svdVals, sample_names = samples, factor_names = as.character(names))
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
countWindowsOverCutoff <- function(qseaSet, GRanges, samples = NULL,
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
#' @param groupMeans Whether to give means of the group column, rather than individual samples.
#' @return A table of data, one row per window
#' @export
#'
getCountTable <- function(qseaSet, groupMeans = FALSE){

  if(groupMeans){

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
#' @param groupMeans Whether to give means of the group column, rather than individual samples.
#' @return A table of data, one row per window
#' @export
#'
getNRPMTable <- function(qseaSet, groupMeans = FALSE){

  if(groupMeans){

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
#' @param groupMeans Whether to give means of the group column, rather than individual samples.
#' @return A table of data, one row per window
#' @export
#'
getBetaTable <- function(qseaSet, groupMeans = FALSE){

  if(groupMeans){

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

#' This function takes a qseaSet and calculates the mean of the nrpm and beta values over the windows given by GRanges.
#' @param qseaSet The qseaSet object.
#' @param GRanges A GRanges object with which windows to take the average of. If missing, all the regions in the qseaSet.
#' @param naMethod What method to use to deal with NA values for the beta values
#' @param minEnrichment The minimum number of reads for a window to be fully methylated, else replace with NA
#' @return A table of data, one row per window
#' @export
#'
getBetaMeans <- function(qseaSet, GRanges = NULL, naMethod = "impute", minEnrichment = 3){

  if(is.null(GRanges)){
    GRanges <- qsea::getRegions(qseaSet)
  }

  samples <- qsea::getSampleNames(qseaSet)

  dataMat <- qseaSet %>%
    filterByOverlaps(GRanges) %>%
    qsea::makeTable(norm_methods = c("nrpm","beta"), samples = samples, minEnrichment = minEnrichment)

  nrpmMat <- dataMat %>%
    dplyr::select(dplyr::matches("nrpm$")) %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "_beta|_nrpm|_means"))

  betaMat <- dataMat %>%
    dplyr::select(dplyr::matches("beta$")) %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "_beta|_nrpm|_means"))

  if(naMethod == "drop"){
    nrpmMat <- nrpmMat %>%
      tidyr::drop_na()

    betaMat <- betaMat %>%
      tidyr::drop_na()
  }

  if(naMethod == "impute"){
    message("Imputing NAs by the mean")
  }

  betaMat %>%
    colMeans(na.rm = TRUE) %>%
    tibble::enframe(name = "sample_name", value = "betaMean") %>%
    dplyr::left_join(nrpmMat %>%
                       colMeans(na.rm = TRUE) %>%
                       tibble::enframe(name = "sample_name", value = "nrpmMean")) %>%
    dplyr::left_join(qsea::getSampleTable(qseaSet)) %>%
    dplyr::mutate(nRegions = nrow(betaMat))

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
    annotateData()

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

setMethod('getSampleNames', 'data.frame',function(object){stop("getSampleNames is not defined on a data frame!")})




#' This function takes a qseaSet and generates a table of data containing data for each sample
#' @param qseaSet The qseaSet object.
#' @param normMethod What normalisation method to use (nrpm or beta)
#' @param groupMeans Number of reads required for a fully methylated region to not be masked (put NA) beta values only.
#' @return A table of data, with a column for each individual sample
#' @export
#'
getDataTable <- function(qseaSet, normMethod = "nrpm", groupMeans = FALSE){

  if(groupMeans){

    qseaSet %>%
      qsea::makeTable(groupMeans =  qsea::getSampleGroups(.), norm_methods = normMethod) %>%
      dplyr::rename_with(~ stringr::str_replace_all(.x, glue::glue("_{normMethod}_means"), "")) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble() %>%
      return()

  } else {

    qseaSet %>%
      qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = normMethod) %>%
      dplyr::rename_with(~ stringr::str_replace_all(.x, glue::glue("_{normMethod}"), "")) %>%
      dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
      tibble::as_tibble() %>%
      return()
  }

}

#' This function takes a qseaSet and writes individual bigWig files with the coverage over each window.
#' @param qseaSet The qseaSet object.
#' @param folderName A folder to save the output bigWig files into
#' @param normMethod What normalisation method to use (e.g. nrpm or beta)
#' @param groupMeans Whether to average over the replicates in the group
#' @param naVal A value to replace NA values with (for beta values).
#' @return A set of bigWig files for each sample in the qseaSet, with coverage over all the windows in the qseaSet
#' @export
#'
writeBigWigs <- function(qseaSet, folderName, normMethod = "nrpm", groupMeans = FALSE, naVal = -1){

  dir.create(folderName, showWarnings = FALSE)

  dataTable <- qseaSet %>%
    getDataTable(normMethod = normMethod, groupMeans = groupMeans)

  if(!groupMeans) {
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

