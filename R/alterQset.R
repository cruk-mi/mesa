#' Filter a qseaSet by overlaps with a GRanges object or table that can be coerced to one
#'
#' This function takes a qseaSet and keeps only a set of regions from it. Either specified as a GRanges object or an index vector of which regions.
#'
#' @param qseaSet The qseaSet object.
#' @param regionsToOverlap A set of windows to keep, either a GRanges object, a dataframe with seqnames, start and end.
#' @return A qseaSet object, with only a subset of the windows.
#' @export
filterByOverlaps <- function(qseaSet, regionsToOverlap){

  if(is.data.frame(regionsToOverlap)) {
    if(length(intersect(colnames(regionsToOverlap),c("seqnames","start","end"))) != 3){
      stop("regionsToOverlap must be a GRanges object or a dataframe with seqnames, start and end.")
    }
  }

  regionsToOverlap <- plyranges::as_granges(regionsToOverlap)

  qseaSetChr <- qseaSet %>% qsea::getRegions() %>% GenomeInfoDb::seqinfo() %>% GenomeInfoDb::seqnames() %>% stringr::str_detect("chr") %>% any()
  windowsChr <- regionsToOverlap %>% GenomeInfoDb::seqinfo() %>% GenomeInfoDb::seqnames() %>% stringr::str_detect("chr") %>% any()

  if(qseaSetChr & !windowsChr){
    regionsToOverlap <- regionsToOverlap %>% tibble::as_tibble() %>% dplyr::mutate(seqnames = paste0("chr",seqnames)) %>% plyranges::as_granges()
  }

  if(!qseaSetChr & windowsChr){
    regionsToOverlap <- regionsToOverlap %>% tibble::as_tibble() %>% dplyr::mutate(seqnames = stringr::str_remove(seqnames,"chr")) %>% plyranges::as_granges()
  }

  qseaSet@count_matrix <- qseaSet@count_matrix[which(plyranges::count_overlaps(qseaSet@regions,regionsToOverlap) > 0), , drop = FALSE]
  qseaSet@regions <- qseaSet@regions %>% plyranges::filter_by_overlaps(regionsToOverlap)

  return(qseaSet)
}

#' Filter a qseaSet by non-overlaps with a GRanges object or table that can be coerced to one
#'
#' This function takes a qseaSet and removes a set of regions from it. Either specified as a GRanges object or a table that can be coerced to one.
#'
#' @param qseaSet The qseaSet object.
#' @param regionsToOverlap A set of windows to keep, either a GRanges object, a dataframe with seqnames, start and end.
#' @return A qseaSet object, with only a subset of the windows.
#' @export
filterByNonOverlaps <- function(qseaSet, regionsToOverlap){

  if(is.data.frame(regionsToOverlap) & ("seqnames" %in% colnames(regionsToOverlap)) &
     ("start" %in% colnames(regionsToOverlap)) & ("end" %in% colnames(regionsToOverlap))){
    stop("regionsToOverlap must be a GRanges object or a dataframe with seqnames, start and end.")
  }

  regionsToOverlap <- plyranges::as_granges(regionsToOverlap)

  qseaSetChr <- qseaSet %>% qsea::getRegions() %>% GenomeInfoDb::seqinfo() %>% GenomeInfoDb::seqnames() %>% stringr::str_detect("chr") %>% any()
  windowsChr <- regionsToOverlap %>% GenomeInfoDb::seqinfo() %>% GenomeInfoDb::seqnames() %>% stringr::str_detect("chr") %>% any()

  if(qseaSetChr & !windowsChr){
    regionsToOverlap <- regionsToOverlap %>% tibble::as_tibble() %>% dplyr::mutate(seqnames = paste0("chr",seqnames)) %>% plyranges::as_granges()
  }

  if(!qseaSetChr & windowsChr){
    regionsToOverlap <- regionsToOverlap %>% tibble::as_tibble() %>% dplyr::mutate(seqnames = stringr::str_remove(seqnames,"chr")) %>% plyranges::as_granges()
  }

  GRangesToKeep <- qsea::getRegions(qseaSet) %>%
    plyranges::filter_by_non_overlaps(plyranges::as_granges(regionsToOverlap))

  qseaSet %>%
    filterByOverlaps(GRangesToKeep) %>%
    return()
}

#' This function takes a qseaSet and adds the information on the number of reads into the sampleTable for easier access.
#'
#' @param qseaSet The qseaSet object.
#' @return A qseaSet object with the sampleTable enhanced with the information on number of reads etc
#' @export
addLibraryInformation <- function(qseaSet){

  curColNames <- qseaSet %>%
    qsea::getSampleTable() %>%
    colnames()

  qseaSet@sampleTable <- cbind(qseaSet@sampleTable,
                               qseaSet@libraries$file_name %>%
                                 dplyr::select(-tidyselect::any_of(curColNames))
  )

  if ("input_file" %in% names(qseaSet@libraries)) {

    curColNames <- qseaSet %>%
      qsea::getSampleTable() %>%
      colnames()

    qseaSet@sampleTable <- cbind(qseaSet@sampleTable,
                                 qseaSet@libraries$input_file %>%
                                   dplyr::rename_with(~ paste0("input_", .x)) %>%
                                   dplyr::select(-tidyselect::any_of(curColNames)) %>%
                                   dplyr::select(-tidyselect::matches("^library_factor$|^offset$"))
    )
  }


  return(qseaSet)
}


#' Subset a qseaSet
#'
#' This function takes a qseaSet and keeps only a subset of the samples from it.
#'
#' @param qseaSet The original qseaSet object.
#' @param samplesToKeep Which samples to keep. Only one of samplesToKeep or samplesToDrop should be specified.
#' @param samplesToDrop Which samples to drop.
#' @return A qseaSet object, with only the selected samples inside it.
#' @export
subsetQset <- function(qseaSet, samplesToKeep = NULL, samplesToDrop = NULL){

  if (length(samplesToKeep) == 0 & length(samplesToDrop) == 0 ) {
   message("No samples remaining, returning an empty qseaSet.")
  }

  if (length(samplesToKeep) > 0 & length(samplesToDrop) > 0 ) {
    stop("Can only specify samples to keep or to drop, not both.")
  }

  samplesNotInSet <- setdiff(c(samplesToKeep,samplesToDrop), qsea::getSampleNames(qseaSet))
  if (length(samplesNotInSet) > 0 ) {
    stop(glue::glue("Sample {samplesNotInSet} not present in the qseaSet!

                    "))
  }

  if (length(samplesToDrop) > 0 ) {
    samplesToKeep = setdiff(qsea::getSampleNames(qseaSet), samplesToDrop)
  }

  newSet <- qseaSet

  newSet@sampleTable <- qseaSet@sampleTable[samplesToKeep,, drop = FALSE]
  newSet@count_matrix <- qseaSet@count_matrix[,samplesToKeep, drop = FALSE]
  newSet@zygosity <- qseaSet@zygosity[samplesToKeep,, drop = FALSE]
  
  if(length(qseaSet@cnv) > 0){
    newSet@cnv <- qseaSet@cnv[,samplesToKeep, drop = FALSE]
  }
  
  newSet@libraries$file_name <- qseaSet@libraries$file_name[samplesToKeep,, drop = FALSE]
  newSet@libraries$input_file <- qseaSet@libraries$input_file[samplesToKeep,, drop = FALSE]
  newSet@enrichment$parameters <- qseaSet@enrichment$parameters[samplesToKeep,, drop = FALSE]
  newSet@enrichment$factors <- qseaSet@enrichment$factors[,samplesToKeep, drop = FALSE]

  return(newSet)
}



#' Rename samples in a qseaSet.
#'
#' This function takes a qseaSet and renames the samples in it.
#'
#' @param qseaSet The qseaSet object.
#' @param pattern Pattern to replace in the names of the samples.
#' @param replacement Pattern to replace with in the names of the samples.
#' @return A qseaSet object, containing all the samples from both qseaSet objects.
#' @export
renameQsetNames <- function(qseaSet, pattern, replacement = "") {

  renamedNames <- qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::pull(sample_name) %>%
    stringr::str_replace_all(pattern, replacement)

  if (all(renamedNames == dplyr::pull(qseaSet@sampleTable, sample_name))) {
    message("Renaming had no effect!")
    return(qseaSet)
  }

  if ( any(duplicated(renamedNames) ) ) {
    stop(glue::glue("Duplicate sample_name now present: {renamedNames[duplicated(renamedNames)]}"))
  }

  newQSet <- qseaSet

  newQSet@sampleTable <- newQSet@sampleTable %>%
    tibble::remove_rownames() %>%
    dplyr::mutate(sample_name = stringr::str_replace_all(sample_name, pattern, replacement),
                  rownameCol = sample_name) %>%
    tibble::column_to_rownames("rownameCol")

  rownames(newQSet@zygosity) <- stringr::str_replace_all(rownames(newQSet@zygosity), pattern, replacement)
  rownames(newQSet@libraries$file_name) <- stringr::str_replace_all(rownames(newQSet@libraries$file_name), pattern, replacement)

  if(!is.null(newQSet@libraries$input_file)) {
    rownames(newQSet@libraries$input_file) <- stringr::str_replace_all(rownames(newQSet@libraries$input_file), pattern, replacement)
  }

  colnames(newQSet@count_matrix) <- stringr::str_replace_all(colnames(newQSet@count_matrix), pattern, replacement)

  rownames(newQSet@enrichment$parameters) <- stringr::str_replace_all(rownames(newQSet@enrichment$parameters), pattern, replacement)
  colnames(newQSet@enrichment$factors) <- stringr::str_replace_all(colnames(newQSet@enrichment$factors), pattern, replacement)

  colnames(GenomicRanges::mcols(newQSet@cnv)) <- stringr::str_replace_all(colnames(GenomicRanges::mcols(newQSet@cnv)), pattern, replacement)

  return(newQSet)
}



#' This function takes a qseaSet and merges samples together into a single sample.
#' @param qseaSet The qseaSet object.
#' @param mergeString A string to merge on
#' @return A qseaSet object with the samples merged together.
#' @export

poolSamples <- function(qseaSet, mergeString){

  ##TODO rewrite this function to use a column.

  samplesToMerge <- qsea::getSampleTable(qseaSet) %>%
    dplyr::pull(sample_name) %>%
    stringr::str_subset(mergeString)

  newNames <- stringr::str_remove(samplesToMerge, mergeString) %>%
    unique()

  if (length(newNames) >= length(samplesToMerge)) {
    stop(glue::glue("New names are the same length as the old ones!"))
  }

  newSet <- qseaSet %>%
    subsetQset(samplesToKeep = samplesToMerge)

  newCounts <- sapply(newNames, function(x) {
    rowSums(qseaSet@count_matrix[,startsWith(colnames(qseaSet@count_matrix), x), drop = FALSE])
  }
  )

  newCNVvals <- sapply(newNames, function(x) {
    #GenomicRanges::mcols(qseaSet@cnv)

    columnsToKeep <- qseaSet %>%
      qsea::getCNV() %>%
      GenomicRanges::mcols() %>%
      colnames() %>%
      stringr::str_subset(x)

    reducedMat <- qseaSet %>%
      qsea::getCNV() %>%
      GenomicRanges::mcols() %>%
      as.matrix()

    weights <- qseaSet@libraries$input_file[columnsToKeep, "total_fragments"]
    weights <- weights/sum(weights)

    # multiply each row by the weight, to average out the number of reads
    rowSums(reducedMat[,columnsToKeep, drop = FALSE] %*% diag(weights) )
  }
  )

  newZygosity <- sapply(newNames, function(x) {

    columnsToKeep <- qseaSet %>%
      qsea::getCNV() %>%
      GenomicRanges::mcols() %>%
      colnames() %>%
      stringr::str_subset(x)

    weights <- qseaSet@libraries$input_file[columnsToKeep, "total_fragments"]
    weights <- weights/sum(weights)

    colSums(qseaSet@zygosity[startsWith(rownames(qseaSet@zygosity), x),,drop = FALSE] * weights)
  }
  ) %>% t()

  newSet@count_matrix <- newCounts
  GenomicRanges::mcols(newSet@cnv) <- newCNVvals
  newSet@zygosity <- newZygosity

  newSampleTable <- purrr::map_dfr(newNames, function(x) {
    qseaSet@sampleTable[startsWith(rownames(qseaSet@sampleTable), x), ,drop = FALSE] %>%
      dplyr::select(-tidyselect::matches("file_name|input_file")) %>%
      dplyr::mutate(
        sample_name = x,
        rownameCol = x) %>%
      dplyr::distinct() %>%
      dplyr::slice(1) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("rownameCol")
  }
  )

  newSet@sampleTable <- newSampleTable

  newLibrariesFile <- purrr::map_dfr(newNames, function(x) {
    qseaSet@libraries$file_name[startsWith(rownames(qseaSet@sampleTable), x),] %>%
      tibble::as_tibble() %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(tidyselect::matches("mean|median|sd|relH|GoGe")), ~stats::weighted.mean(.x,total_fragments)),
                       dplyr::across(tidyselect::any_of(tidyselect::matches("reads|pairs|r1s|fragments")), ~sum(.x))
      ) %>%
      dplyr::mutate(rownameCol = x, library_factor = NA, offset = NA) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("rownameCol")

  }
  )

  newLibrariesInput <- purrr::map_dfr(newNames, function(x) {
    qseaSet@libraries$input_file[startsWith(rownames(qseaSet@sampleTable), x),] %>%
      tibble::as_tibble() %>%
      dplyr::mutate(total_fragments2 = total_fragments) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(tidyselect::matches("mean|median|sd|relH|GoGe")), ~stats::weighted.mean(.x,total_fragments)),
                       dplyr::across(tidyselect::any_of(tidyselect::matches("reads|pairs|r1s|fragments")), ~sum(.x))
      ) %>%
      dplyr::mutate(rownameCol = x, library_factor = NA, offset = NA) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("rownameCol")

  }
  )

  newSet@libraries$file_name <- as.matrix(newLibrariesFile)
  newSet@libraries$input_file <- as.matrix(newLibrariesInput)

  oldSet <- qseaSet %>%
    subsetQset(samplesToDrop = samplesToMerge)

  finalSet <-  combineQsets(oldSet, newSet) %>%
    addNormalisation(enrichmentMethod = "blind1-15")

  return(finalSet)

}

#' This function takes a qseaSet and renames the samples based on a column of the sampleTable.
#' @param qseaSet The qseaSet object.
#' @param newNameColumn A column of the qseaSet to use to rename each sample with.
#' @return A qseaSet object with the samples renamed.
#' @export

renameSamples <- function(qseaSet, newNameColumn){

  newNameColumn <- rlang::enquo(newNameColumn)

  renamedNames <- qseaSet %>% qsea::getSampleTable() %>% dplyr::pull(!!newNameColumn)

  if (any(is.na(renamedNames))) {
    stop(glue::glue("NAs present in the new sample name column"))
  }

  if (all(renamedNames == dplyr::pull(qseaSet@sampleTable, sample_name))) {
    message("Renaming had no effect!")
    return(qseaSet)
  }
  if (any(duplicated(renamedNames))) {
    stop(glue::glue("Duplicate sample_name now present: {renamedNames[duplicated(renamedNames)]}"))
  }
  newQSet <- qseaSet
  newQSet@sampleTable <- newQSet@sampleTable %>% tibble::remove_rownames() %>%
    dplyr::mutate(sample_name = renamedNames,
                  rownameCol = renamedNames) %>%
    tibble::column_to_rownames("rownameCol")

  rownames(newQSet@zygosity) <- renamedNames
  rownames(newQSet@libraries$file_name) <- renamedNames
  if (!is.null(newQSet@libraries$input_file)) {
    rownames(newQSet@libraries$input_file) <- renamedNames
  }
  colnames(newQSet@count_matrix) <- renamedNames
  rownames(newQSet@enrichment$parameters) <- renamedNames
  colnames(newQSet@enrichment$factors) <- renamedNames
  colnames(GenomicRanges::mcols(newQSet@cnv)) <- renamedNames
  return(newQSet)
}
