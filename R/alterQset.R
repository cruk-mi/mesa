#' Subset a qseaSet by overlaps with genomic regions
#'
#' Retain only windows that overlap a set of genomic regions. This is a light
#' wrapper around [plyranges::filter_by_overlaps()] that also harmonizes
#' chromosome naming (with/without `"chr"` prefix) between inputs.
#'
#' @param qseaSet A `qseaSet` object.
#' @param regionsToOverlap Genomic regions to use for filtering. Either a
#'   [GenomicRanges::GRanges()] or a `data.frame` with columns
#'   `seqnames`, `start`, and `end` (coerced internally to `GRanges`).
#'
#' @return A filtered `qseaSet` containing only windows that overlap
#'   `regionsToOverlap`.
#'
#' @details
#' If one input uses `chr1` while the other uses `1`, the function will adjust
#' the region seqnames internally so that overlaps are computed correctly.
#'
#' @seealso
#'   [plyranges::filter_by_overlaps()], [filterByNonOverlaps()],
#'   [GenomicRanges::GRanges()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Using a data.frame
#' rgns <- data.frame(seqnames = 7, start = 25_002_001, end = 25_017_900)
#' qs_sub <- filterByOverlaps(exampleTumourNormal, regionsToOverlap = rgns)
#' qs_sub
#'
#' # Using a GRanges
#' gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(25_002_001, 25_017_900))
#' qs_sub2 <- filterByOverlaps(exampleTumourNormal, regionsToOverlap = gr)
#' qs_sub2
#'
#' @rdname alterQsetOverlap
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


#' Subset a qseaSet by non‑overlaps with genomic regions
#'
#' Retain only windows that **do not** overlap a set of genomic regions. This is
#' a wrapper around [plyranges::filter_by_non_overlaps()] with automatic
#' chromosome naming harmonization.
#'
#' @return A filtered `qseaSet` containing only windows that do **not** overlap
#'   `regionsToOverlap`.
#'
#' @seealso
#'   [plyranges::filter_by_non_overlaps()], [filterByOverlaps()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' rgns <- data.frame(seqnames = 7, start = 25_002_001, end = 25_017_900)
#' qs_sub <- filterByNonOverlaps(exampleTumourNormal, regionsToOverlap = rgns)
#' qs_sub
#'
#' @rdname alterQsetOverlap
#' @export
filterByNonOverlaps <- function(qseaSet, regionsToOverlap){

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

  GRangesToKeep <- qsea::getRegions(qseaSet) %>%
    plyranges::filter_by_non_overlaps(plyranges::as_granges(regionsToOverlap))

  qseaSet %>%
    filterByOverlaps(GRangesToKeep) %>%
    return()
}


#' Add library metrics into the qseaSet sample table
#'
#' Convenience helper to append columns from `qseaSet@libraries$file_name`
#' (and, when available, from `qseaSet@libraries$input_file` prefixed with
#' `"input_"`) to the `sampleTable`. Existing columns are not duplicated.
#'
#' @param qseaSet A `qseaSet` object.
#'
#' @return The same `qseaSet`, with an augmented `sampleTable`.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- addLibraryInformation(exampleTumourNormal)
#' head(qsea::getSampleTable(qs))
#'
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


#' Subset a qseaSet by samples
#'
#' Keep only a specified set of samples (or, alternatively, drop a set). All
#' internal matrices/slots are subset consistently.
#'
#' @param qseaSet A `qseaSet` object.
#' @param samplesToKeep Character vector of sample names to keep. Use **either**
#'   `samplesToKeep` **or** `samplesToDrop`.
#' @param samplesToDrop Character vector of sample names to drop. Use **either**
#'   `samplesToKeep` **or** `samplesToDrop`.
#'
#' @return A `qseaSet` containing only the selected samples.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' subsetQset(exampleTumourNormal, samplesToKeep = c("Colon1_T","Colon1_N"))
#'
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


#' Rename samples in a qseaSet by regex
#'
#' Apply a regex replacement to sample names and update all relevant slots
#' consistently (counts, CNV, enrichment, zygosity, libraries, sampleTable).
#'
#' @param qseaSet A `qseaSet` object.
#' @param pattern A regex `pattern` to match in sample names.
#' @param replacement A replacement string (default: `""`).
#'
#' @return A `qseaSet` with renamed samples.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' # Replace "T" with "Tumour" in sample names
#' qs2 <- renameQsetNames(exampleTumourNormal, pattern = "T$", replacement = "Tumour")
#' head(colnames(qs2@count_matrix))
#'
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


#' Pool (merge) samples that share a common name prefix
#'
#' Merge samples whose names share a common prefix by removing a suffix pattern
#' (`mergeString`) and summing counts across the resulting groups. CNV and
#' zygosity are combined using fragment‐count weights. A new `qseaSet` is
#' returned with pooled samples and updated metadata.
#'
#' @param qseaSet A `qseaSet` object.
#' @param mergeString A regex used to remove the **suffix** that distinguishes
#'   samples to be merged. The remaining prefix becomes the pooled sample name.
#'   After removal, pooled names are matched using `startsWith()`.  
#'   For example, with names like `Patient1_T` and `Patient1_N`, use
#'   `mergeString = "_[TN]$"` to pool into `Patient1`.
#'
#' @return A `qseaSet` with pooled samples.
#'
#' @details
#' Internally, pooled counts are simple sums. CNV values and zygosity are
#' combined by a weighted average using `total_fragments` from
#' `@libraries$input_file`. Library statistics are aggregated (sums for counts,
#' weighted means for rates).
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#' # Pool Tumour/Normal pairs by patient into a single sample per patient:
#' # e.g., "Colon1_T","Colon1_N" -> "Colon1"
#' qs_pool <- poolSamples(exampleTumourNormal, mergeString = "_[TN]$")
#' qs_pool
#' }
#'
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

  newSet@libraries$file_name <- as.data.frame(newLibrariesFile)
  newSet@libraries$input_file <- as.data.frame(newLibrariesInput)

  oldSet <- qseaSet %>%
    subsetQset(samplesToDrop = samplesToMerge)

  finalSet <-  combineQsets(oldSet, newSet) %>%
    addNormalisation(enrichmentMethod = "blind1-15")

  return(finalSet)

}


#' Rename samples using a column from the sample table
#'
#' Replace sample names with values from a user‑supplied column in
#' `qseaSet@sampleTable`, updating all dependent slots consistently.
#'
#' @param qseaSet A `qseaSet` object.
#' @param newNameColumn Name (unquoted) of a column in the sample table that
#'   contains the new sample names. Must be unique and contain no `NA`s.
#'
#' @return A `qseaSet` with sample names replaced by `newNameColumn`.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal %>%
#'   dplyr::mutate(newName = paste0("Sample", seq_len(nrow(qsea::getSampleTable(exampleTumourNormal)))))
#' qs_renamed <- renameSamples(qs, newNameColumn = newName)
#' head(qsea::getSampleTable(qs_renamed)$sample_name)
#'
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
