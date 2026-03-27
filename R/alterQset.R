#' Subset a qseaSet by overlaps / non-overlaps with genomic regions
#'
#' `filterByOverlaps()` retains windows that **overlap**, while
#' `filterByNonOverlaps()` retains windows that **do not overlap** a set of
#' genomic regions. Inputs may use different chromosome styles (with or without
#' the `"chr"` prefix); styles are harmonized automatically before overlap
#' computation.
#'
#' @param qseaSet `qseaSet`
#'   A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param regionsToOverlap `GRanges` or `data.frame`
#'   Genomic regions used for filtering. If a `data.frame`, it must contain
#'   columns `seqnames` (`character()`), `start` (`integer()`), and `end`
#'   (`integer()`); the object is coerced internally to a `GRanges`.
#'
#' @return A `qseaSet` object:
#'
#' * **filterByOverlaps()**: returns the input `qseaSet` restricted to
#'   windows that **overlap** `regionsToOverlap`.
#'
#' * **filterByNonOverlaps()**: returns the input `qseaSet` restricted to
#'   windows that **do not overlap** `regionsToOverlap`.
#'
#' @details
#' Chromosome naming is aligned using [`GenomeInfoDb::seqlevelsStyle()`] so that
#' overlap computation is consistent even when one input uses `"chr1"` and the
#' other uses `"1"`. Only the **window set** is filtered; sample metadata and
#' counts are subset accordingly.
#'
#' @seealso 
#' \code{\link[plyranges]{filter_by_overlaps}},
#' \code{\link[plyranges]{filter_by_non_overlaps}},
#' \code{\link[GenomicRanges]{GRanges-class}}
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Data frame input
#' rgns <- data.frame(seqnames = 7L, start = 25002001L, end = 25017900L)
#' filterByOverlaps(exampleTumourNormal, regionsToOverlap = rgns)
#' filterByNonOverlaps(exampleTumourNormal, regionsToOverlap = rgns)
#'
#' # GRanges input
#' gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(25002001L, 25017900L))
#' filterByOverlaps(exampleTumourNormal, regionsToOverlap = gr)
#' filterByNonOverlaps(exampleTumourNormal, regionsToOverlap = gr)
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


#' Add library metrics to the qseaSet sample table
#'
#' Append library information to the `sampleTable` of a `qseaSet`. Columns are
#' taken from `qseaSet@libraries$file_name`, and when available,
#' from `qseaSet@libraries$input_file` (prefixed with `"input_"`). Existing
#' columns are not duplicated.
#'
#' @param qseaSet `qseaSet`
#'   A qseaSet object containing sequencing data.
#'
#' @return A `qseaSet` object:
#'
#' * **addLibraryInformation()**: returns the same `qseaSet` with its
#'   `sampleTable` augmented by library-level columns.
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
#' Restrict a `qseaSet` to a specified set of samples. All internal
#' matrices and slots are subset consistently.  
#' Use exactly one of `samplesToKeep` or `samplesToDrop`.
#'
#' @param qseaSet `qseaSet`  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param samplesToKeep `character()`  
#'   Names of samples to retain. Cannot be used together with `samplesToDrop`.
#'   **Default:** `NULL`.
#'   
#' @param samplesToDrop `character()`  
#'   Names of samples to remove. Cannot be used together with `samplesToKeep`.
#'   **Default:** `NULL`.
#'   
#' @return A `qseaSet` object:
#'
#' * **subsetQset()**: returns the input `qseaSet` restricted to the
#'   selected samples.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Keep only two samples
#' subsetQset(exampleTumourNormal, samplesToKeep = c("Colon1_T", "Colon1_N"))
#'
#' # Drop one sample
#' subsetQset(exampleTumourNormal, samplesToDrop = "Lung1_N")
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
    samplesToKeep <- setdiff(qsea::getSampleNames(qseaSet), samplesToDrop)
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
#' Apply a regular expression replacement to sample names and update all
#' relevant slots consistently (counts, CNV, enrichment, zygosity, libraries,
#' and the `sampleTable`).
#'
#' @param qseaSet `qseaSet`  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param pattern `character(1)`  
#'   Regular expression used to match substrings in sample names.
#'
#' @param replacement `character(1)`  
#'   String to replace matches of `pattern`. Defaults to `""` (empty string).
#'
#' @return A `qseaSet` object:
#'
#' * **renameQsetNames()**: returns the input `qseaSet` with sample names
#'   updated according to `pattern` and `replacement`. All dependent slots
#'   are updated consistently.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Replace "T" with "Tumour" in sample names
#' renameQsetNames(exampleTumourNormal, 
#'                 pattern = "T",
#'                 replacement = "Tumour")
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

  asValidNames <- base::make.names(renamedNames)
  
  if(any(asValidNames != renamedNames)){
    stop(glue::glue("Sample names must be valid names for columns in R without quoting.
  See the help for base::make.names, but generally use only letters, numbers, 
  underscores and dots, and names can't start with a number. 
  Issues were found with: 
    {paste(renamedNames[renamedNames != asValidNames], collapse = '\n    ')}
   "))
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
#' zygosity values are combined using fragment‐count weights. The result is a
#' new `qseaSet` with pooled samples and updated metadata.
#'
#' @param qseaSet `qseaSet`  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param mergeString `character(1)`  
#'   Regular expression used to remove the **suffix** that distinguishes samples
#'   to be merged. The remaining prefix becomes the pooled sample name. After
#'   removal, pooled names are matched using `startsWith()`.  
#'   For example, with names like `"Patient1_T"` and `"Patient1_N"`, use  
#'   `mergeString = "_[TN]$"` to pool them into `"Patient1"`.
#'
#' @return A `qseaSet` object:
#'
#' * **poolSamples()**: returns a new `qseaSet` where samples with names sharing
#'   a common prefix (after applying `mergeString`) are merged. Counts are summed,
#'   CNV and zygosity values are combined using weighted averages, and library
#'   statistics are aggregated.
#'
#' @details
#' Internally, pooled counts are simple sums. CNV values and zygosity are
#' combined by a weighted average using `total_fragments` from
#' `@libraries$input_file`. Library statistics are aggregated (sums for counts,
#' weighted means for rates).
#'
#' @examples
#' \dontrun{
#' data(exampleTumourNormal, package = "mesa")
#' # Pool Tumour/Normal pairs by patient into a single sample per patient:
#' # e.g., "Colon1_T","Colon1_N" -> "Colon1"
#' poolSamples(exampleTumourNormal, mergeString = "_[TN]$")
#' }
#' @export
poolSamples <- function(qseaSet, mergeString){

  ##TODO: rewrite this function to use a column.
  ##TODO: check this actually works! 
  ##TODO: make sure it works in general

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

  newCounts <- vapply(
    newNames,
    function(x) {
      rowSums(
        qseaSet@count_matrix[
          , startsWith(colnames(qseaSet@count_matrix), x),
          drop = FALSE
        ]
      )
    },
    FUN.VALUE = numeric(nrow(qseaSet@count_matrix))
  )

  newCNVvals <- vapply(
    newNames,
    function(x) {
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
      weights <- weights / sum(weights)
      
      rowSums(
        reducedMat[, columnsToKeep, drop = FALSE] %*% diag(weights)
      )
    },
    FUN.VALUE = numeric(nrow(GenomicRanges::mcols(qseaSet@cnv)))
  )

  newZygosity <- vapply(
    newNames,
    function(x) {
      columnsToKeep <- qseaSet %>%
        qsea::getCNV() %>%
        GenomicRanges::mcols() %>%
        colnames() %>%
        stringr::str_subset(x)
      
      weights <- qseaSet@libraries$input_file[columnsToKeep, "total_fragments"]
      weights <- weights / sum(weights)
      
      colSums(
        qseaSet@zygosity[
          startsWith(rownames(qseaSet@zygosity), x),
          ,
          drop = FALSE
        ] * weights
      )
    },
    FUN.VALUE = numeric(ncol(qseaSet@zygosity))
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
#' Replace sample names with values from a user-supplied column in the
#' `sampleTable` of a `qseaSet`, updating all dependent slots consistently.
#'
#' @param qseaSet `qseaSet`  
#'   A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param newNameColumn `character(1)`  
#'   The name of a column in the `sampleTable` containing new sample names.
#'   Values must be unique and contain no `NA`s.  
#'   (Unquoted column names are supported via tidy evaluation.)
#'
#' @return A `qseaSet` object:
#'
#' * **renameSamples()**: returns the input `qseaSet` with sample names replaced
#'   by values from `newNameColumn`. All associated slots are updated to reflect
#'   the new naming.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' 
#' exampleTumourNormal %>%
#'   dplyr::mutate(newName = paste0("Sample", seq_len(nrow(qsea::getSampleTable(.))))) %>%
#'   renameSamples(newNameColumn = "newName")
#'   
#' @export
renameSamples <- function(qseaSet, newNameColumn){

  newNameColumn <- rlang::enquo(newNameColumn)

  renamedNames <- qseaSet %>% qsea::getSampleTable() %>% dplyr::pull(!!newNameColumn)

  if (any(is.na(renamedNames))) {
    stop(glue::glue("NAs present in the new sample name column"))
  }

  asValidNames <- base::make.names(renamedNames)
  
  if(any(asValidNames != renamedNames)){
    stop(glue::glue("Sample names must be valid names for columns in R without quoting.
  See the help for base::make.names, but generally use only letters, numbers, 
  underscores and dots, and names can't start with a number. 
  Issues were found with: 
    {paste(renamedNames[renamedNames != asValidNames], collapse = '\n    ')}
   "))
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
