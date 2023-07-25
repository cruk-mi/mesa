#' Fit a generalised linear model
#'
#' This function fits a generalised linear model given a qseaSet, and then add a set of contrasts to it
#'
#' @param qseaSet The qseaSet object
#' @param variable The variable to treat as the independent variable between samples in the linear model
#' @param covariates A vector including any additional covariates to include in the model
#' @param formula A formula to use
#' @param contrasts A data frame with the comparisons to do
#' @param keepIndex A vector of indices to keep, overwrites minReadCount.
#' @param minReadCount A minimum read count for a row to be considered, in the qseaSet as provided.
#' @param minNRPM A minimum normalised reads per million to apply
#' @param checkPVals Whether to check excessive numbers of the p-values are exactly zero, to catch a bug in qsea.
#' @param calcDispersionAll Whether to use samples that are not present in the contrasts to fit the initial generalised linear model, including them in the calculation of dispersion estimates.
#' Setting this to be TRUE will mean that adding additional samples to the qseaSet will change the calculated DMRs, even if they are not being compared across.
#' @return A qseaGLM object
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
    numExtraSamples <- length(setdiff(getSampleNames(qseaSet, samplesInContrasts)))
    if(numExtraSamples > 0){
        message(glue::glue("Calculating dispersion estimates including {numExtraSamples} that are not being used in contrasts."))
    }
  }

  if (is.null(keepIndex) & minNRPM == 0 ) {
    keepIndex = which(matrixStats::rowMaxs(qsea::getCounts(qseaSet)) >= minReadCount)
  }

  if (is.null(keepIndex) & minNRPM >= 0) {

    keepIndex <- qseaSet %>%
      qsea::makeTable(samples = samplesInContrasts,
                      norm_methods = "nrpm",
                      chunksize = 1000000, verbose = FALSE) %>%
      dplyr::select(tidyselect::matches("nrpm")) %>%
      apply(1,max) %>%
      {which(. >= minNRPM)}

    names(keepIndex) <- NULL

  }

  # make a design object based on the formula
  design <- stats::model.matrix(formula, qseaSet %>% qsea::getSampleTable())


  message(glue::glue("Fitting initial GLM on {length(keepIndex)} windows, using {BiocParallel::bpworkers()} cores"))

  qseaGLM <- suppressMessages(qsea::fitNBglm(qseaSet,
                                             design,
                                             keep = keepIndex,
                                             minRowSum = 0,
                                             norm_method = "beta",
                                             parallel = (BiocParallel::bpworkers() > 1),
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
                                                  parallel = TRUE,
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


#' Extract the data table from all the contrasts
#'
#' This function extracts the information from the DMR contrasts into a wide table
#'
#' @param qseaSet The qseaSet object
#' @param qseaGLM A qseaGLM object
#' @param sampleNames Which samples to return the information for.
#' @param variable Which variable the DMR comparison has been performed on
#' @param fdrThres False discovery rate threshold.
#' @param keepData Whether to keep the individual data columns
#' @param keepGroupMeans Whether to keep the group means
#' @param keepPvals Whether to keep the unadjusted p-values in the output
#' @param keepFragmentInfo Whether to keep information on the average fragment length and MAPQ (if present)
#' @param direction Whether to use regions that are up/down/both.
#' @return A tibble with the data
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

  if (length(sigIndex) == 0) { print("No windows found.") }

  #TODO Quick hack to stop error when a contrast has one or no significantly varying windows
  hack <- FALSE
  if (length(sigIndex) == 0 ) {
    print("No windows found, using hack to return, needs testing")
    sigIndex <- 1
    hack <- TRUE
  }

  if (length(sigIndex) == 1 ) {
    print("One window found, using hack to return, needs testing")
    sigIndex <- c(sigIndex,sigIndex)
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

  # Remove the added window from the hack
  if (hack) {
    dataTable <- dataTable[-1,]
  }

  return(dataTable)
}

#' Make a full set of possible contrasts for calculating DMRs on.
#'
#' This function returns a set of all possible contrasts
#'
#' @param qseaSet The qseaSet object
#' @param variable Which variable to use to calculate the DMRs between
#' @return A tibble with two columns, each row with a pair of contrasts
#' @export
#'
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

#' Fit the GLM and return the data
#'
#' This function fits a GLM and returns the resulting data table
#'
#' @param qseaSet The qseaSet object
#' @param variable Which variable to use to calculate the DMRs between
#' @param covariates Any variables to use as covariates, for instance patient in a paired analysis
#' @param contrasts A data frame with two columns, group1 and group2, with the strings to compare between in each. Multiple rows means that multiple comparisons will be fitted
#' @param formula Alternative formula mode for calculating DMRs (not recommended)
#' @param minNRPM A minimum normalised reads per million value that at least one sample (in the contrasts) must reach in order to consider the region for further calculation. Set this or minReadCount (but preferably this).
#' @param minReadCount A minimum read count that at least one sample (in the contrasts) must reach in order to consider the region for further calculation. Preferably use minNRPM.
#' @param fdrThres False discovery rate threshold.
#' @param keepContrastMeans Whether to keep the columns containing the means of the contrasts in the output
#' @param keepData Whether to keep the individual data columns in the output
#' @param keepGroupMeans Whether to keep the group means in the output
#' @param keepPvals Whether to keep the unadjusted p-values in the output
#' @param checkPVals Whether to check that the p-values aren't mostly zero to avoid a bug with covariates, only turn this off if you are sure what you are doing!
#' @param direction Whether to keep regions that are up/down/both.#'
#' @param calcDispersionAll Whether to use samples that are not present in the contrasts to fit the initial generalised linear model, including them in the calculation of dispersion estimates.
#' Setting this to be TRUE will mean that adding additional samples to the qseaSet will change the calculated DMRs, even if they are not being compared across.
#' @return A tibble with the data
#' @export
#'
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
