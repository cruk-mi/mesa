#' Mix two samples to generate a synthetic qseaSet sample
#'
#' Create a new synthetic sample in a `qseaSet` by mixing reads from two
#' existing samples in user-specified proportions. This can be useful for
#' benchmarking, downsampling experiments, or simulating mixtures of tumour
#' and normal samples.
#'
#' @param qseaSet A `qseaSet` object.
#' @param sample1 `character(1)` First sample name, contributing
#'   `proportion * nReadsTotal` reads.
#' @param sample2 `character(1)` Second sample name, contributing the remainder
#'   of reads.
#' @param nReadsTotal `integer(1)` Total number of reads to simulate in the
#'   new mixed sample.
#' @param proportion `numeric(1)` Proportion of reads taken from `sample1`.
#'   Must be between 0 and 1. The remaining reads are taken from `sample2`.
#' @param newName `character(1)` Name for the new synthetic sample. If `NULL`,
#'   a name is generated automatically.
#' @param groupName `character(1)` Group name for the new sample in the
#'   `sampleTable`. Defaults to `newName`.
#' @param onlyNew `logical(1)` If `TRUE`, return only the synthetic sample.
#'   If `FALSE` (default), return the original `qseaSet` combined with the new
#'   sample.
#' @param renormalise `logical(1)` Whether to run [addNormalisation()] on the
#'   result. For efficiency, set to `FALSE` if repeatedly mixing, and run
#'   normalisation once at the end.
#'
#' @return A `qseaSet` containing the new synthetic sample. By default, the
#'   original samples are preserved as well.
#'
#' @seealso [mixThreeQsetSamples()], [downSample()]
#' @family sample-simulation
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#'
#' ## Mix two tumour replicates 70/30, keeping original samples
#' qs_mix <- mixSamples(qs, "LUAD_1", "LUAD_2", nReadsTotal = 1e5,
#'                      proportion = 0.7, newName = "LUAD_mix")
#'
#' ## Return only the synthetic mixture sample
#' qs_mix_only <- mixSamples(qs, "LUAD_1", "LUAD_2", nReadsTotal = 5e4,
#'                           proportion = 0.5, onlyNew = TRUE)
#'
#' @export
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


#' Mix three samples to generate a synthetic qseaSet sample
#'
#' Create a new synthetic sample in a `qseaSet` by mixing reads from three
#' existing samples according to user-specified proportions. This function
#' is currently intended for internal use and has not been extensively tested.
#'
#' @param qseaSet A `qseaSet` object.
#' @param sample1, sample2, sample3 `character(1)` Sample names to mix.
#' @param nReadsTotal `integer(1)` Total number of reads in the new sample.
#' @param proportion1 `numeric(1)` Proportion of reads taken from `sample1`.
#'   Must be between 0 and 1.
#' @param proportion2 `numeric(1)` Proportion of reads taken from `sample2`.
#'   Must be between 0 and 1.
#'
#'   The proportion of reads from `sample3` is computed as
#'   `1 - proportion1 - proportion2`. The sum of `proportion1 + proportion2`
#'   must be strictly less than 1.
#'
#' @param newName `character(1)` Name for the new synthetic sample. If `NULL`,
#'   a name is generated automatically.
#' @param groupName `character(1)` Group name for the new sample in the
#'   `sampleTable`. Defaults to `newName`.
#' @param onlyNew `logical(1)` If `TRUE`, return only the synthetic sample.
#'   If `FALSE` (default), return the original `qseaSet` combined with the new
#'   sample.
#' @param renormalise `logical(1)` Whether to run [addNormalisation()] on the
#'   result. For efficiency, set to `FALSE` if repeatedly mixing, and run
#'   normalisation once at the end.
#'
#' @return A `qseaSet` containing the new synthetic sample.
#'
#' @seealso [mixSamples()]
#' @family sample-simulation
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#'
#' ## Mix three samples: 50% tumour, 30% normal, 20% third sample
#' qs_mix3 <- mixThreeQsetSamples(qs, "LUAD_1", "Normal_1", "LUAD_2",
#'                                nReadsTotal = 1e5,
#'                                proportion1 = 0.5, proportion2 = 0.3,
#'                                newName = "ThreeMix")
#'
#' ## Only return the synthetic mixture
#' qs_mix3_only <- mixThreeQsetSamples(qs, "LUAD_1", "Normal_1", "LUAD_2",
#'                                     nReadsTotal = 5e4,
#'                                     proportion1 = 0.6, proportion2 = 0.2,
#'                                     onlyNew = TRUE)
#'
#' @keywords internal
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


