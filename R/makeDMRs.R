#' Fit a generalized linear model (GLM) to a qseaSet
#'
#' Fit a negative-binomial GLM using [qsea::fitNBglm()] on counts from a
#' [qsea::qseaSet], and add one or more contrasts.  
#' This is the low-level function used by [calculateDMRs()].
#'
#' @param qseaSet A [qsea::qseaSet] object containing counts and sample metadata.
#' @param variable Character scalar. Column in the sample table used as the
#'   primary explanatory variable.
#' @param covariates Character vector of additional covariates (e.g. batch,
#'   patient ID).
#' @param formula Optional `formula` overriding `variable`/`covariates`.
#' @param contrasts Data frame with columns `group1` and `group2`, specifying
#'   contrasts to compute.
#' @param keepIndex Optional integer vector of row indices to retain; overrides
#'   filtering by `minReadCount`/`minNRPM`.
#' @param minReadCount Minimum raw count required in at least one sample to
#'   retain a window.
#' @param minNRPM Minimum NRPM required in at least one sample to retain a window.
#' @param checkPVals Logical. If `TRUE`, stop or warn if excessive proportions
#'   of p-values are exactly zero (indicating possible instability).
#' @param calcDispersionAll Logical. Whether to use samples that are not present
#' in the contrasts to fit the initial generalised linear model, including them 
#' in the calculation of dispersion estimates. Setting this to be TRUE will mean
#' that adding additional samples to the qseaSet will change the calculated DMRs,
#' even if they are not being compared across.
#' @return A `qseaGLM` object with fitted model and contrasts.
#'
#' @details
#' - Dispersion is estimated across all contrasts (or all samples, if
#'   `calcDispersionAll = TRUE`).  
#' - Contrasts are added sequentially to the fitted GLM using
#'   [qsea::addContrast()].  
#' - By default, beta-normalised counts are used as the outcome.
#'
#' @family DMR-detection
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' contr <- tibble::tibble(group1 = "LUAD", group2 = "NormalLung")
#' glmfit <- fitQseaGLM(qs, variable = "type", contrasts = contr)
#' glmfit
#' 
#' @export
fitQseaGLM <- function(qseaSet, variable = NULL,  covariates = NULL,
                       contrasts = NULL, keepIndex = NULL, minReadCount = 0, minNRPM = 1,
                       checkPVals = TRUE, formula = NULL, calcDispersionAll = FALSE){

  if (!is.null(contrasts)) {
    nContrasts <- nrow(contrasts)
  }else{
    nContrasts <- 0
  }

  if (nContrasts > 0 & is.null(variable) & is.data.frame(contrasts)) {

    if ("variable" %in% colnames(contrasts)) {

      variablesInTable <- contrasts %>%
        dplyr::pull(variable) %>%
        unique()

      if (length(variablesInTable) == 1) {
        variable <- variablesInTable
      } else {
        stop(glue::glue("Multiple variables present in contrasts: {variablesInTable} "))
      }

    } else if (!is.null(formula)) {
      variable <- all.vars(formula)[1]
    } else {
      stop(glue::glue("Do not know what variable to use! Provide in the contrast list or as an argument."))
    }
  }

  # Make a progress bar
  pb <- progress::progress_bar$new(total = (nContrasts + 1))

  #message("Fitting full GLM w/ adjustments")

  if (is.null(formula)) {
    formula <- stats::as.formula(paste0("~", paste(c(variable, covariates), collapse = "+"), " + 0"))
  }

  if (!all(covariates %in% colnames( qsea::getSampleTable(qseaSet)))) {
    stop(glue::glue("Covariate {setdiff(covariates,colnames( qsea::getSampleTable(qseaSet)))} missing from the sampleTable."))
  }

  contrasts <- contrasts %>%
    dplyr::select(-tidyselect::matches("^variable$"))

  if ( ncol(contrasts) == 2) {
    colnames(contrasts) <- c("group1", "group2")
  } else {
    colnames(contrasts)[1:2] <- c("group1", "group2")
  }

  valuesInContrasts <- contrasts %>%
    tidyr::pivot_longer(tidyselect::starts_with("group"), names_to = "name", values_to = "group") %>%
    dplyr::pull(group)

  samplesInContrasts <-  qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::filter(!!rlang::sym(variable) %in% valuesInContrasts) %>%
    dplyr::pull(sample_name)

  if(!calcDispersionAll){
    qseaSet <- qseaSet %>%
      filter(sample_name %in% samplesInContrasts)

  } else {
    numExtraSamples <- length(setdiff(qsea::getSampleNames(qseaSet), samplesInContrasts))
    if(numExtraSamples > 0){
        message(glue::glue("Calculating dispersion estimates including {numExtraSamples} samples that are not being used in contrasts."))
    }
  }

  if (is.null(keepIndex) & minNRPM == 0 ) {
    keepIndex = which(matrixStats::rowMaxs(qsea::getCounts(qseaSet %>% dplyr::filter(sample_name %in% samplesInContrasts))) >= minReadCount)
  }

  if (is.null(keepIndex) & minNRPM >= 0) {

    keepIndex <- qseaSet %>%
      getDataTable(normMethod = "nrpm",
                   addMethodSuffix = TRUE) %>%
      dplyr::select(tidyselect::matches("nrpm")) %>%
      apply(1,max) %>%
      {which(. >= minNRPM)}

    names(keepIndex) <- NULL

  }

    # make a design object based on the formula
  design <- stats::model.matrix(formula, qseaSet %>% qsea::getSampleTable())

  if(getMesaParallel()){
    message(glue::glue("Fitting initial GLM on {length(keepIndex)} windows, using {BiocParallel::bpworkers()} cores"))
  } else {
    message(glue::glue("Fitting initial GLM on {length(keepIndex)} windows, without using parallelisation."))
  }

  qseaGLM <- suppressMessages(qsea::fitNBglm(qseaSet,
                                             design,
                                             keep = keepIndex,
                                             minRowSum = 0,
                                             norm_method = "beta",
                                             parallel = getMesaParallel(),
                                             verbose = FALSE))

  pb$tick()

  # Yes, a for loop. The issue is that it adds repeatedly to the qseaGLM object, so can't be vectorised easily.
  for (i in seq_along(1:nrow(contrasts))) {

    conName <- paste0(variable, contrasts[i,"group1"], "-", variable, contrasts[i,"group2"])

    if (!(contrasts[i,"group1"] %in% qsea::getSampleTable(qseaSet)[,variable])) {
      stop(glue::glue("value {contrasts[i,]$group1} not found in column {variable} of the sampleTable!"))
    }

    if (!(contrasts[i,"group2"] %in% qsea::getSampleTable(qseaSet)[,variable])) {
      stop(glue::glue("value {contrasts[i,]$group2} not found in column {variable} of the sampleTable!"))
    }

    if ("name" %in% colnames(contrasts)) {
      conNameClean <- contrasts[i,"name"] %>% dplyr::pull()
    }else{
      # Remove the hyphen from the name, because it messes up things later.
      # Also remove the variable name (whatever it is)
      conNameClean <- conName %>%
        stringr::str_replace("-","_vs_") %>%
        stringr::str_remove_all(variable)
    }

    message(glue::glue("Performing contrast {conNameClean}"))

    limContrast <- limma::makeContrasts(contrasts = conName, levels = design)

    qseaGLM <- suppressMessages(qsea::addContrast(qseaSet,
                                                  qseaGLM,
                                                  contrast = limContrast,
                                                  name = conNameClean,
                                                  parallel = getMesaParallel(),
                                                  verbose = FALSE))

    pb$tick()
    if (mean(qseaGLM@contrast[[conNameClean]]$LRT_pval == 0) >= 0.2 & checkPVals) {

      if(is.null(covariates)){
        warning("Warning! More than 20% of windows have p-values of exactly 0, possibly something has gone wrong! \n
           Set checkPVals = FALSE to ignore this.")
      } else {

        stop("Error! More than 20% of windows have p-values of exactly 0, possibly an error! \n
           Try not including covariates in the model if included, or set checkPVals = FALSE to ignore this if sure.")
      }
    }
  }
  message("Contrasts Complete")
  return(qseaGLM)
}


#' Extract DMR-level results from a fitted GLM
#'
#' Given a [qseaGLM] object with one or more contrasts, extract a wide table of
#' statistics, group means, and (optionally) per-sample values.
#'
#' @param qseaSet A [qsea::qseaSet] object (used for metadata and counts).
#' @param qseaGLM A fitted [qseaGLM] object (from [fitQseaGLM()]).
#' @param sampleNames Character vector of sample names to include; if `NULL`,
#'   determined by `keepData`.
#' @param variable Variable on which contrasts were based (for group means).
#' @param fdrThres Numeric(1). FDR threshold for selecting significant regions.
#' @param keepData Logical. If `TRUE`, include per-sample data columns.
#' @param keepGroupMeans Logical. If `TRUE`, include group mean columns.
#' @param keepPvals Logical. If `TRUE`, include raw (unadjusted) p-values.
#' @param keepFragmentInfo Logical. If `TRUE`, include fragment length / MAPQ info.
#' @param direction Character: `"up"`, `"down"`, or `"both"`. Passed to
#'   [qsea::isSignificant()].
#'
#' @return A tibble with one row per DMR and columns for statistics,
#'   adjusted p-values, and (optionally) group means or per-sample values.
#'
#' @family DMR-detection
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' contr <- tibble::tibble(group1 = "LUAD", group2 = "NormalLung")
#' glmfit <- fitQseaGLM(qs, variable = "type", contrasts = contr)
#' dmrTab <- getDMRsData(qs, glmfit, variable = "type", fdrThres = 0.1)
#' head(dmrTab)
#' 
#' @export
getDMRsData <- function(qseaSet, qseaGLM, sampleNames = NULL, variable = NULL, keepData = FALSE, keepGroupMeans = FALSE,
                        fdrThres = 0.05, keepPvals = FALSE, keepFragmentInfo = FALSE,
                        direction = "both"){

  sampleTable <- qsea::getSampleTable(qseaSet)

  if (is.null(sampleNames)){
    if(keepData) {
      sampleNames <- sampleTable$sample_name
    } else {
      sampleNames <- NULL
    }
  }

  sigIndex <- purrr::map(names(qseaGLM@contrast), ~ qsea::isSignificant(qseaGLM, contrast = ., fdr_th = fdrThres, direction = direction)) %>%
    unlist() %>%
    unique()

  hack <- FALSE
  # workaround to ensure that qsea doesn't drop the data frame to a vector when it subsets, copied row removed at the end
  if (length(sigIndex) == 1 ) {
    sigIndex <- c(sigIndex,sigIndex)
    hack <- TRUE
  }

  # Use the same workaround to ensure that the same columns and types are returned as would be expected if any DMRs were found.
  if (length(sigIndex) == 0 ) {
    sigIndex <- 1
    hack <- TRUE
  }

  if (!is.null(variable)) {

    contrastMeansList <- sampleTable %>%
      dplyr::pull(!!variable) %>%
      unique() %>%
      stats::na.omit() %>%
      stats::setNames(nm = .) %>%
      purrr::map(function(x) {sampleTable %>%
          dplyr::filter(!!rlang::sym(variable) == x) %>%
          dplyr::pull(sample_name)})

  }

  if (keepGroupMeans & variable != "group") {

    groupMeansList <- getSampleGroups2(qseaSet)[names(getSampleGroups2(qseaSet)) %in% unique(sampleTable$group)] %>%
      c(contrastMeansList)

  } else {
    groupMeansList <- contrastMeansList
  }

  dataTable <- qsea::makeTable(qseaSet, qseaGLM,
                               keep = sigIndex,
                               norm_methods = c("beta","nrpm"),
                               samples = if (keepData) {sampleNames} else{NULL},
                               groupMeans = groupMeansList, verbose = FALSE)

  if (!keepPvals) {
    dataTable <- dataTable %>%
      dplyr::select(-tidyselect::matches("pvalue"))
  }

  if (!keepFragmentInfo) {
    dataTable <- dataTable %>%
      dplyr::select(-tidyselect::matches("avgFragment"))
  }

  # Remove the added window from the hacky workaround
  if (hack) {
    dataTable <- dataTable[-1,]
  }

  return(dataTable)
}


#' Generate all possible pairwise contrasts
#'
#' Construct a tibble of all pairwise contrasts between levels of a categorical
#' variable in the sample table of a [qsea::qseaSet].
#'
#' @param qseaSet A [qsea::qseaSet] object.
#' @param variable Character scalar. Column in the sample table.
#'
#' @return A tibble with columns `group1` and `group2`.
#'
#' @family DMR-detection
#' 
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' makeAllContrasts(exampleTumourNormal, "type")
#' 
#' @export
makeAllContrasts <- function(qseaSet, variable){
  vals <- qseaSet %>%
    qsea::getSampleTable() %>%
    tidyr::drop_na(tidyselect::all_of(variable)) %>%
    dplyr::pull(variable) %>%
    unique() %>%
    gtools::mixedsort()

  tidyr::expand_grid(group1 = vals, group2 = vals) %>%
    filter(group1 < group2) %>%
    return()
}


#' Fit GLM and return DMR data in one step
#'
#' High-level wrapper that calls [fitQseaGLM()] and [getDMRsData()] to produce a
#' DMR results table directly.
#'
#' @inheritParams fitQseaGLM
#' @inheritParams getDMRsData
#' @param contrasts Data frame or string specifying contrasts. If `"All"`,
#'   generate all pairwise contrasts via [makeAllContrasts()]. Strings of the
#'   form `"A_vs_B"`, `"All_vs_X"`, or `"X_vs_All"` are also supported.
#' @param keepContrastMeans Logical. If `FALSE`, remove contrast mean columns
#'   from the output.
#'
#' @return A tibble with one row per significant region, containing DMR
#'   statistics, adjusted p-values, and delta-beta values for each contrast.
#'
#' @family DMR-detection
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' # One contrast
#' calculateDMRs(qs, variable = "type", contrasts = "LUAD_vs_NormalLung")
#'
#' # All pairwise contrasts
#' dmrTab <- calculateDMRs(qs, variable = "type", contrasts = "All")
#' head(dmrTab)
#' 
#' @export
calculateDMRs <- function(qseaSet,
                          variable = NULL,
                          covariates = NULL,
                          contrasts = NULL,
                          minReadCount = 0,
                          minNRPM = 1,
                          checkPVals = TRUE,
                          fdrThres = 0.05,
                          keepPvals = FALSE,
                          formula = NULL,
                          keepContrastMeans = TRUE,
                          keepData = FALSE,
                          keepGroupMeans = FALSE,
                          direction = "both",
                          calcDispersionAll = FALSE){

  if (is.null(variable)) {stop("variable must be specified!")}
  if (is.null(contrasts)) {stop("contrasts must be specified!")}

  if (!is.data.frame(contrasts)) {
    if (contrasts %in% c("All", "all")) {
      contrasts <- makeAllContrasts(qseaSet, variable)
      message(
        glue::glue(
          "Calculating {nrow(contrasts)} possible contrasts on the {variable} column."
        )
      )
    } else if (contrasts %in% c("First", "first")) {
      contrasts <- makeAllContrasts(qseaSet, variable)[1, ]
      message(glue::glue(
        "Calculating the first possible contrast on the {variable} column."
      ))
    } else if (stringr::str_detect(contrasts,"All_vs_|all_vs_")){

      value2 = contrasts %>% stringr::str_remove("All_vs_|all_vs_")
      contrasts <- tibble::tibble(group1 = qseaSet %>% pull(variable) %>% unique() %>% setdiff(value2),
                                  group2 = value2)
      message(glue::glue(
        "Calculating all ({nrow(contrasts)}) possible contrasts against {value2} on the {variable} column."
      ))
    } else if (stringr::str_detect(contrasts,"_vs_All")){

      value1 = contrasts %>% stringr::str_remove("_vs_All|_vs_all")
      contrasts <- tibble::tibble(group1 = value1,
                                  group2 = qseaSet %>% pull(variable) %>% unique() %>% setdiff(value1))
      message(glue::glue(
        "Calculating all ({nrow(contrasts)}) possible contrasts between {value1} and the rest of {variable} column."
      ))
    } else if (stringr::str_detect(contrasts,"_vs_")){
      value1 <- stringr::str_remove(contrasts,"_vs_.*")
      value2 <- stringr::str_remove(contrasts,".*_vs_")
      contrasts <- tibble::tibble(group1 = value1, group2 = value2)
    } else {
      stop(glue::glue("String {contrasts} not recognised."))
    }

  }

  contrasts <- contrasts %>%
    dplyr::select(-tidyselect::matches("^variable$"))

  if( ncol(contrasts) == 2){
    colnames(contrasts) <- c("group1", "group2")
  } else {
    stop("Contrasts data frame should contain columns group1 and group2 (or exactly two columns).")
  }

  if(is.null(contrasts)){stop("No contrasts specified!")}

  qseaGLM <- fitQseaGLM(qseaSet, variable = variable,  covariates = covariates,
                        contrasts = contrasts,  minReadCount = minReadCount,
                        minNRPM = minNRPM,
                        checkPVals = checkPVals, formula = formula,
                        calcDispersionAll = calcDispersionAll)

  dataTable <- getDMRsData(qseaSet, qseaGLM, sampleNames = qsea::getSampleNames(qseaSet),
                           fdrThres = fdrThres, keepPvals = keepPvals, keepData = keepData,
                           keepGroupMeans = keepGroupMeans,
                           variable = variable,
                           direction = direction) %>%
    dplyr::rename(seqnames = chr, start = window_start, end = window_end)

  deltas  <- purrr::map_dfc(1:nrow(contrasts),
                            function(x){
                              name1 <- contrasts[x,]$group1
                              name2 <- contrasts[x,]$group2
                              (dataTable[,paste0(name1,"_beta_means")] - dataTable[,paste0(name2,"_beta_means")]) %>%
                                tibble::enframe(name = "rowIndex", value = paste0(name1,"_vs_",name2,"_deltaBeta")) %>%
                                dplyr::select(-rowIndex)
                            }
  )

  dataTable <- dataTable %>%
    dplyr::bind_cols(deltas) %>%
    tibble::as_tibble()

  # ewww, a for loop. Moves the deltaBeta columns around.
  for(adjPvalString in (dataTable %>% colnames() %>% stringr::str_subset("_adjPval$"))){
    dataTable <- dataTable %>% dplyr::relocate(stringr::str_replace(adjPvalString, "_adjPval$","_deltaBeta"), .after = !!adjPvalString)
  }

  if (!keepContrastMeans) {

    contrastNames <- contrasts %>% {c(pull(.,group1), pull(.,group2))}
    colsToRemove <- paste0(contrastNames, rep(c("_beta_means","_nrpm_means"),rep(length(contrastNames),2)))

    dataTable <- dataTable %>%
      dplyr::select(-tidyselect::all_of(colsToRemove))
  }

  return(dataTable)

}
