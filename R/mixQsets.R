#' Mix two samples to generate a synthetic qseaSet sample
#'
#' Create a new synthetic sample in a `qseaSet` by mixing reads from two
#' existing samples in user-specified proportions. Useful for benchmarking,
#' downsampling experiments, or simulating tumour–normal mixtures.
#' 
#' Note that if the qseaSet windows have been filtered prior to using this, 
#' the fraction of fragments present in the current qseaSet is calculated, 
#' such that the total number of reads would have been approximately `nReadsTotal`.
#' E.g. if there are 50% of the total fragments in the windows present in the current samples, 
#' 0.5*`nReadsTotal` fragments will be sampled.
#' 
#'
#' @param qseaSet `qseaSet`.  
#'   The input object containing the two source samples.
#'
#' @param sample1 `character(1)`.  
#'   First sample name, contributing `proportion * nReadsTotal` reads.
#'
#' @param sample2 `character(1)`.  
#'   Second sample name, contributing the remainder of reads.
#'
#' @param nReadsTotal `integer(1)`.  
#'   Total number of reads to simulate in the new mixed sample.  
#'   Must be `>= 0`.
#'
#' @param proportion `numeric(1)`.  
#'   Proportion of reads taken from `sample1`. Must be in `[0, 1]`.  
#'   The remaining `1 - proportion` reads are taken from `sample2`.
#'
#' @param newName `character(1)` or `NULL`.  
#'   Name for the new synthetic sample.  
#'   **Default:** `NULL` (a name is generated automatically, e.g.,
#'   `"Mix_<sample1>_<sample2>_<proportion>"`).
#'
#' @param groupName `character(1)` or `NULL`.  
#'   Group name for the new sample recorded in the `sampleTable`.  
#'   **Default:** `NULL` (uses `newName`).
#'
#' @param onlyNew `logical(1)`.  
#'   If `TRUE`, return only the synthetic sample. If `FALSE`, return the
#'   original `qseaSet` with the new sample appended.  
#'   **Default:** `FALSE`.
#'
#' @param renormalise `logical(1)`.  
#'   Whether to run [addNormalisation()] on the result. For efficiency, set
#'   to `FALSE` when creating multiple mixtures and run normalisation once at
#'   the end.  
#'   **Default:** `TRUE`.
#'
#' @return A `qseaSet` containing the new synthetic sample:  
#' * If `onlyNew = TRUE`, a `qseaSet` with **one** (mixed) sample.  
#' * If `onlyNew = FALSE`, the input `qseaSet` **plus** the mixed sample.  
#' The `sampleTable`, counts, and relevant slots are updated consistently.
#'
#' @details
#' * Input validation ensures `sample1`/`sample2` exist, `proportion ∈ [0,1]`,
#'   and `nReadsTotal` is non-negative.  
#' * The mixing strategy (e.g., proportional allocation vs sampling) follows
#'   the package’s implementation; set `renormalise = TRUE` to recompute
#'   offsets/enrichment after adding the new sample if required by your workflow.
#'
#' @seealso
#' [mixThreeQsetSamples()], [downSample()], [addNormalisation()]
#'
#' @family sample-simulation
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' 
#' 
#' # Mix two samples 50/50 into ~100k total fragments; append to the qseaSet
#' exampleTumourNormal %>% 
#'   mixSamples("Colon1_T", 
#'              "Colon1_N", 
#'              nReadsTotal = 100000, 
#'              proportion = 0.5, 
#'              renormalise = FALSE) %>%
#'                            
#' # Only a few reads are included in this small subset of windows (82):
#'  getCountTable() %>% 
#'  pull(Mix_Colon1_T_Colon1_N_0.5) %>% 
#'  sum()
#'
#'
#' # Mix two samples 80/20 and return only the synthetic sample, with an explicit name and group
#' exampleTumourNormal %>%
#'   mixSamples("Colon1_T", "Colon1_N",
#'              nReadsTotal = 50000,
#'              proportion  = 0.8,
#'              newName     = "Mix_Colon1_T_0.8_Colon1_N_0.2",
#'              groupName   = "Synthetic",
#'              onlyNew     = TRUE ,
#'              renormalise = FALSE) %>%
#'              
#'   # Only a few reads are included in this small subset of windows (45):
#'   getCountTable() %>% 
#'   dplyr::pull(Mix_Colon1_T_0.8_Colon1_N_0.2) %>% 
#'   sum()
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

  windows1 <- sample(rep(seq_along(nrow(counts)), counts[,sample1]), replace = FALSE, size = nReads1 * onTargetFrac1)
  windows2 <- sample(rep(seq_along(nrow(counts)), counts[,sample2]), replace = FALSE, size = nReads2 * onTargetFrac2)

  newCounts <- c(windows1, windows2) %>%
    tibble::enframe(name = NULL, value = "window") %>%
    dplyr::count(window, name = "n") %>%
    dplyr::left_join(tibble::tibble(window = seq_along(nrow(counts))),., copy = TRUE, by = "window") %>%
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
#' Create a new synthetic sample in a `qseaSet` by mixing reads from **three**
#' existing samples according to user-specified proportions. Useful for simple
#' mixture simulations; currently intended for internal use.
#' 
#' Note that if the qseaSet windows have been filtered prior to using this,
#' the fraction of fragments present in the current qseaSet is calculated, 
#' such that the total number of reads would have been approximately `nReadsTotal`.
#' E.g. if there are 50% of the total fragments in the windows present in the current samples, 
#' 0.5*`nReadsTotal` fragments will be sampled.
#'
#' @param qseaSet `qseaSet`.  
#'   The input object containing the three source samples.
#'
#' @param sample1, sample2, sample3 `character(1)`.  
#'   Sample names to mix (must exist in `qseaSet`).
#'
#' @param nReadsTotal `integer(1)`.  
#'   Total number of reads to simulate in the mixed sample. Must be `>= 0`.
#'
#' @param proportion1 `numeric(1)`.  
#'   Proportion taken from `sample1`. Must be in `[0, 1]`.  
#'   **Default:** none (must be supplied).
#'
#' @param proportion2 `numeric(1)`.  
#'   Proportion taken from `sample2`. Must be in `[0, 1]`.  
#'   **Default:** none (must be supplied).
#'
#' The proportion from `sample3` is computed as `1 - proportion1 - proportion2`.  
#' The sum `proportion1 + proportion2` must be **strictly less than 1** (so that
#' `sample3` gets a non-negative share).
#'
#' @param newName `character(1)` or `NULL`.  
#'   Name for the new synthetic sample.  
#'   **Default:** `NULL` (auto-generated, e.g., `"Mix3_<s1>_<s2>_<s3>_<p1>_<p2>"`).
#'
#' @param groupName `character(1)` or `NULL`.  
#'   Group name recorded in the `sampleTable` for the new sample.  
#'   **Default:** `NULL` (uses `newName`).
#'
#' @param onlyNew `logical(1)`.  
#'   If `TRUE`, return only the synthetic mixture. If `FALSE`, append it to the
#'   original `qseaSet`.  
#'   **Default:** `FALSE`.
#'
#' @param renormalise `logical(1)`.  
#'   Whether to run [addNormalisation()] on the result. For batching multiple
#'   mixtures, consider `FALSE` and normalise once at the end.  
#'   **Default:** `TRUE`.
#'
#' @return A `qseaSet` containing the new synthetic sample:  
#' * If `onlyNew = TRUE`, a `qseaSet` with **one** (mixed) sample.  
#' * If `onlyNew = FALSE`, the input `qseaSet` **plus** the mixed sample.  
#' The `sampleTable`, counts, and dependent slots are updated consistently.
#'
#' @details
#' * Validates that `sample1`, `sample2`, `sample3` exist; `nReadsTotal >= 0`;
#'   `proportion1, proportion2 ∈ [0, 1]`; and `proportion1 + proportion2 ≤ 1`.  
#' * The mixing strategy follows the package implementation (proportional
#'   allocation / sampling); use `renormalise = TRUE` if downstream steps rely
#'   on offsets/enrichment being recomputed.
#'
#' @seealso
#' [mixSamples()] (two-sample mixture), [addNormalisation()]
#'
#' @family sample-simulation
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' 
#' 
#' # Mix three samples 50/30/20 into ~100k total fragments; append to the qseaSet
#' exampleTumourNormal %>% 
#'   mixThreeQsetSamples("Colon1_T", "Colon1_N", "Lung1_T",
#'                        nReadsTotal = 100000,
#'                        proportion1 = 0.5, proportion2 = 0.3,
#'                        newName = "Mix_Colon1_T_0.5_Colon1_N_0.3_Lung1_T_0.2",
#'                        renormalise = FALSE) %>%
#'                        
#'  # Only a few reads are included in this small subset of windows (78):
#'  getCountTable() %>% 
#'  pull(Mix_Colon1_T_0.5_Colon1_N_0.3_Lung1_T_0.2) %>% 
#'  sum()
#'
#'
#' # Mix three samples and return only the synthetic sample, with an explicit name and group.
#' exampleTumourNormal %>%
#'   mixThreeQsetSamples("Colon1_T", "Colon1_N", "Lung1_T",
#'              nReadsTotal = 50000,
#'              proportion1 = 0.1, proportion2 = 0.2,
#'              newName     = "Mix_Colon1_T_0.1_Colon1_N_0.2_Lung1_T_0.7",
#'              groupName   = "Synthetic",
#'              onlyNew     = TRUE ,
#'              renormalise = FALSE) %>%
#'              
#'   # Only a few reads are included in this small subset of windows (29):
#'   getCountTable() %>% 
#'   dplyr::pull(Mix_Colon1_T_0.1_Colon1_N_0.2_Lung1_T_0.7) %>% 
#'   sum()
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

  windows1 <- sample(rep(seq_along(nrow(counts)), counts[,sample1]), replace = FALSE, size = nReads1 * onTargetFrac1)
  windows2 <- sample(rep(seq_along(nrow(counts)), counts[,sample2]), replace = FALSE, size = nReads2 * onTargetFrac2)
  windows3 <- sample(rep(seq_along(nrow(counts)), counts[,sample3]), replace = FALSE, size = nReads3 * onTargetFrac3)

  newCounts <- c(windows1, windows2, windows3) %>%
    table() %>%
    tibble::enframe(name = "window") %>%
    dplyr::mutate(window = as.integer(window)) %>%
    dplyr::left_join(tibble::tibble(window = seq_along(nrow(counts))),., copy = TRUE, by = "window") %>%
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
