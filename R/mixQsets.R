#' This function takes a qseaSet and makes a new sample by mixing two samples
#' @param qseaSet The qseaSet object.
#' @param sample1 First sample name, from which to take proportion of samples
#' @param sample2 Second sample name
#' @param nReadsTotal Number of reads in total to have after mixing
#' @param proportion The proportion to take from sample1, the rest will come from sample2
#' @param newName A name to give the new sample
#' @param groupName A name to use in the group column in the sampleTable
#' @param onlyNew Whether to only return the new sample.
#' @param renormalise Whether to renormalise the result. Speeds up the process when you are repeatedly subsampling, only need to do it once at the end.
#' @return A qseaSet object with an extra
#' @export
#'
mixSamples <- function(qseaSet, sample1, sample2, nReadsTotal, proportion, newName = NULL, groupName = NULL,
                           onlyNew = FALSE,
                           renormalise = TRUE){

  qsea:::checkSamples(qseaSet,c(sample1, sample2))

  if(is.null(newName)) { newName <- paste0("Mix","_",sample1,"_",sample2,"_",proportion) }

  message(newName)

  if(is.null(groupName)) { groupName <- newName }

  nReads1 <- ceiling(nReadsTotal * proportion)
  nReads2 <- nReadsTotal - nReads1

  if (proportion < 0 || proportion > 1) {
    stop(glue::glue("Proportion must be a number between 0 and 1, not {proportion}."))
  }

  validFragNums <- qseaSet@libraries$file_name[c(sample1, sample2),"valid_fragments", drop = FALSE] %>%
    tibble::rownames_to_column() %>%
    tibble::deframe()

  counts <- qseaSet@count_matrix[,c(sample1, sample2)]
  sampleColSums <- colSums(counts)

  message(glue::glue("Mixing {nReads1} reads from {sample1} with {nReads2} reads from {sample2}."))

  onTargetFrac1 <- sampleColSums[sample1]/validFragNums[sample1]
  onTargetFrac2 <- sampleColSums[sample2]/validFragNums[sample2]

  if (sampleColSums[sample1]  < nReads1 * onTargetFrac1) {
    stop(glue::glue("Not enough reads in {sample1}, only {sampleColSums[sample1]} out of {nReads1 * onTargetFrac1} requested."))
  }

  if (sampleColSums[sample2]  < nReads2 * onTargetFrac2) {
    stop(glue::glue("Not enough reads in {sample2}, only {sampleColSums[sample2]} out of {nReads2 * onTargetFrac2} requested."))
  }

  windows1 <- sample(rep(1:nrow(counts), counts[,sample1]), replace = FALSE, size = nReads1 * onTargetFrac1)
  windows2 <- sample(rep(1:nrow(counts), counts[,sample2]), replace = FALSE, size = nReads2 * onTargetFrac2)

  newCounts <- c(windows1, windows2) %>%
    tibble::enframe(name = NULL, value = "window") %>%
    dplyr::count(window, name = "n") %>%
    dplyr::left_join(tibble::tibble(window = 1:nrow(counts)),., copy = TRUE, by = "window") %>%
    dplyr::mutate(n = tidyr::replace_na(n,0)) %>%
    dplyr::select(n)

  newSet <- qseaSet %>%
    subsetQset(samplesToKeep = qsea::getSampleNames(qseaSet)[1]) %>%
    renameQsetNames(paste0("^",qsea::getSampleNames(qseaSet)[1],"$"), newName)

  newSet@count_matrix <- as.matrix(newCounts)
  colnames(newSet@count_matrix) <- newName

  reducedCNVMat <- qseaSet %>%
    qsea::getCNV() %>%
    GenomicRanges::mcols() %>%
    as.matrix()

  GenomicRanges::mcols(newSet@cnv) <- proportion * reducedCNVMat[,sample1] + (1 - proportion) * reducedCNVMat[,sample2]
  colnames(GenomicRanges::mcols(newSet@cnv)) <- newName

  newSet@zygosity <- proportion * qseaSet@zygosity[sample1,,drop = FALSE] + (1 - proportion) * qseaSet@zygosity[sample2,,drop = FALSE]
  rownames(newSet@zygosity) <- newName

  newSet@sampleTable <- tibble::tibble(sample_name = newName,
                                       rownameCol = newName,
                                       group = groupName,
                                       type = qseaSet %>% qsea::getSampleTable() %>% dplyr::filter(sample_name == !!sample1) %>% dplyr::pull(type),
                                       tumour = qseaSet %>% qsea::getSampleTable() %>% dplyr::filter(sample_name == !!sample1) %>% dplyr::pull(tumour),
                                       sample1 = sample1,
                                       sample2 = sample2,
                                       prop1 = proportion) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("rownameCol") %>%
    as.data.frame()

  newSet@libraries$file_name[,"valid_fragments"] <- nReadsTotal
  newSet@libraries$file_name[,"offset"] <- NA
  newSet@libraries$file_name[,"library_factor"] <- NA

  rownames(newSet@libraries$file_name) <- newName

  if(renormalise){
    newSet <- addNormalisation(newSet, enrichmentMethod = "blind1-15")
  }

  if (onlyNew) {
    return(newSet)
  } else {
    return(combineQsets(qseaSet, newSet))
  }

}



#' This function takes a qseaSet and makes a new sample by mixing three samples. Currently internal only as untested.
#' @param qseaSet The qseaSet object.
#' @param sample1 First sample name, from which to take proportion of samples
#' @param sample2 Second sample name
#' @param sample3 Third sample name
#' @param nReadsTotal Number of reads in total to have after mixing
#' @param proportion1 The proportion to take from sample1
#' @param proportion2 The proportion to take from sample2
#' @param newName A name to give the new sample
#' @param groupName A name to use in the group column in the sampleTable
#' @param onlyNew Whether to only return the new sample.
#' @param renormalise Whether to renormalise the result. Speeds up the process when you are repeatedly subsampling, only need to do it once at the end.
#' @return A qseaSet object with an extra
#'
mixThreeQsetSamples <- function(qseaSet, sample1, sample2, sample3, nReadsTotal, proportion1, proportion2, newName = NULL, groupName = NULL,
                                onlyNew = FALSE,
                                renormalise = TRUE){

  qsea:::checkSamples(qseaSet,c(sample1, sample2, sample3))

  if(is.null(newName)) { newName <- paste0("Mix","_",sample1,"_",sample2,"_",sample3,"_",proportion1,"_",proportion2) }

  message(newName)

  if(is.null(groupName)) { groupName <- newName }

  nReads1 <- ceiling(nReadsTotal * proportion1)
  nReads2 <- ceiling(nReadsTotal * proportion2)
  nReads3 <- nReadsTotal - nReads1 - nReads2

  if (proportion1 < 0 || proportion1 > 1) {
    stop(glue::glue("Proportion1 must be a number between 0 and 1, not {proportion1}."))
  }

  if (proportion2 < 0 || proportion2 > 1) {
    stop(glue::glue("Proportion2 must be a number between 0 and 1, not {proportion2}."))
  }


  if (proportion1 + proportion2 >= 1) {
    stop(glue::glue("Sum of proportions must be less than 1, not {proportion1 + proportion2}."))
  }

  validFragNums <- qseaSet@libraries$file_name[c(sample1, sample2, sample3),"valid_fragments", drop = FALSE] %>%
    tibble::rownames_to_column() %>%
    tibble::deframe()

  counts <- qseaSet@count_matrix[,c(sample1, sample2, sample3)]
  sampleColSums <- colSums(counts)

  message(glue::glue("Mixing {nReads1} reads from {sample1} with {nReads2} reads from {sample2} and {nReads3} reads from {sample3}."))

  nReads1 <- ceiling(nReadsTotal * proportion1)
  nReads2 <- ceiling(nReadsTotal * proportion2)
  nReads3 <- nReadsTotal - nReads1 - nReads2

  onTargetFrac1 <- sampleColSums[sample1]/validFragNums[sample1]
  onTargetFrac2 <- sampleColSums[sample2]/validFragNums[sample2]
  onTargetFrac3 <- sampleColSums[sample3]/validFragNums[sample3]

  if (sampleColSums[sample1]  < nReads1 * onTargetFrac1) {
    stop(glue::glue("Not enough reads in {sample1}, only {sampleColSums[sample1]} out of {nReads1 * onTargetFrac1} requested."))
  }

  if (sampleColSums[sample2]  < nReads2 * onTargetFrac2) {
    stop(glue::glue("Not enough reads in {sample2}, only {sampleColSums[sample2]} out of {nReads2 * onTargetFrac2} requested."))
  }

  if (sampleColSums[sample3]  < nReads3 * onTargetFrac2) {
    stop(glue::glue("Not enough reads in {sample3}, only {sampleColSums[sample3]} out of {nReads3 * onTargetFrac3} requested."))
  }

  windows1 <- sample(rep(1:nrow(counts), counts[,sample1]), replace = FALSE, size = nReads1 * onTargetFrac1)
  windows2 <- sample(rep(1:nrow(counts), counts[,sample2]), replace = FALSE, size = nReads2 * onTargetFrac2)
  windows3 <- sample(rep(1:nrow(counts), counts[,sample3]), replace = FALSE, size = nReads3 * onTargetFrac3)

  newCounts <- c(windows1, windows2, windows3) %>%
    table() %>%
    tibble::enframe(name = "window") %>%
    dplyr::mutate(window = as.integer(window)) %>%
    dplyr::left_join(tibble::tibble(window = 1:nrow(counts)),., copy = TRUE, by = "window") %>%
    dplyr::mutate(value = tidyr::replace_na(value,0)) %>%
    dplyr::select(value)

  newSet <- qseaSet %>%
    subsetQset(samplesToKeep = qsea::getSampleNames(qseaSet)[1]) %>%
    renameQsetNames(paste0("^",qsea::getSampleNames(qseaSet)[1],"$"), newName)

  newSet@count_matrix <- as.matrix(newCounts)
  colnames(newSet@count_matrix) <- newName

  reducedCNVMat <- qseaSet %>%
    qsea::getCNV() %>%
    GenomicRanges::mcols() %>%
    as.matrix()

  GenomicRanges::mcols(newSet@cnv) <- proportion1 * reducedCNVMat[,sample1] + proportion2 * reducedCNVMat[,sample2] + (1 - proportion1 - proportion2) * reducedCNVMat[,sample3]
  colnames(GenomicRanges::mcols(newSet@cnv)) <- newName

  newSet@zygosity <- proportion1 * qseaSet@zygosity[sample1,,drop = FALSE] + proportion2 * qseaSet@zygosity[sample2,,drop = FALSE] + (1 - proportion1 - proportion2) * qseaSet@zygosity[sample3,,drop = FALSE]
  rownames(newSet@zygosity) <- newName

  newSet@sampleTable <- tibble::tibble(sample_name = newName,
                                       rownameCol = newName,
                                       group = groupName,
                                       type = qseaSet %>% qsea::getSampleTable() %>% dplyr::filter(sample_name == !!sample1) %>% dplyr::pull(type),
                                       tumour = qseaSet %>% qsea::getSampleTable() %>% dplyr::filter(sample_name == !!sample1) %>% dplyr::pull(tumour),
                                       sample1 = sample1,
                                       sample2 = sample2,
                                       prop1 = proportion1,
                                       prop2 = proportion2) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("rownameCol") %>%
    as.data.frame()

  newSet@libraries$file_name[,"valid_fragments"] <- nReadsTotal
  newSet@libraries$file_name[,"offset"] <- NA
  newSet@libraries$file_name[,"library_factor"] <- NA

  rownames(newSet@libraries$file_name) <- newName

  if(renormalise){
    newSet <- addNormalisation(newSet, enrichmentMethod = "blind1-15")
  }

  if (onlyNew) {
    return(newSet)
  } else {
    return(combineQsets(qseaSet, newSet))
  }

}


