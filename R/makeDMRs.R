#' Fit a generalised linear model
#'
#' This function fits a generalised linear model given a qseaSet, and then add a set of contrasts to it
#'
#' @param qseaSet The qseaSet object
#' @param variable The variable to treat as the independent variable between samples in the linear model
#' @param covariates A vector including any additional covariates to include in the model
#' @param formula A formula to use
#' @param contrastsToDo A data frame with the comparisons to do
#' @param keepIndex A vector of indices to keep, overwrites minReadCount.
#' @param minReadCount A minimum read count for a row to be considered, in the qseaSet as provided.
#' @param minNRPM A minimum normalised reads per million to apply
#' @param checkPVals Whether to check excessive numbers of the p-values are exactly zero, to catch a bug in qsea.
#' @return A qseaGLM object
#' @export
fitQseaGLM <- function(qseaSet, variable = NULL,  covariates = NULL,
                       contrastsToDo = NULL, keepIndex = NULL, minReadCount = 0, minNRPM = 1,
                       checkPVals = TRUE, formula = NULL){

  sampleTable <- qsea::getSampleTable(qseaSet)

  if (!is.null(contrastsToDo)) {
    nContrasts <- nrow(contrastsToDo)
  }else{
    nContrasts <- 0
  }

  if (nContrasts > 0 & is.null(variable) & is.data.frame(contrastsToDo)) {

    if ("variable" %in% colnames(contrastsToDo)) {

      variablesInTable <- contrastsToDo %>%
        dplyr::pull(variable) %>%
        unique()

      if (length(variablesInTable) == 1) {
        variable <- variablesInTable
      } else {
        stop(glue::glue("Multiple variables present in contrastsToDo: {variablesInTable} "))
      }

    } else if (!is.null(formula)) {
      variable <- all.vars(formula)[1]
    } else {
      stop(glue::glue("Do not know what variable to use! Provide in the contrast list or as an argument."))
    }
  }

  # Make a progress bar
  pb <- progress::progress_bar$new(total = (nContrasts + 1))

  message("Fitting full GLM w/ adjustments")

  if (is.null(formula)) {
    formula <- stats::as.formula(paste0("~", paste(c(variable, covariates), collapse = "+"), " + 0"))
  }

  if (!all(covariates %in% colnames(sampleTable))) {
    stop(glue::glue("Covariate {setdiff(covariates,colnames(sampleTable))} missing from the sampleTable."))
  }

  # make a design object based on the formula
  design <- stats::model.matrix(formula, sampleTable)

  valuesInContrasts <- contrastsToDo %>%
    tidyr::pivot_longer(tidyselect::starts_with("sample"), names_to = "name", values_to = "sample") %>%
    dplyr::pull(sample)

  samplesInContrasts <- sampleTable %>%
    dplyr::filter(!!rlang::sym(variable) %in% valuesInContrasts) %>%
    dplyr::pull(sample_name)

  if (is.null(keepIndex) & minNRPM == 0 ) {
    keepIndex = which(matrixStats::rowMaxs(qsea::getCounts(qseaSet)[,samplesInContrasts]) >= minReadCount)
  }

  if (is.null(keepIndex) & minNRPM >= 0) {

    keepIndex <- qseaSet %>%
      qsea::makeTable(samples = samplesInContrasts,
                      norm_methods = "nrpm", chunksize = 1000000) %>%
      dplyr::select(tidyselect::matches("nrpm")) %>%
      apply(1,max) %>%
      {which(. >= minNRPM)}

    names(keepIndex) <- NULL

  }

  print("Fitting initial GLM")
  qseaGLM <- suppressMessages(qsea::fitNBglm(qseaSet,
                                             design,
                                             keep = keepIndex,
                                             minRowSum = 0,
                                             norm_method = "beta",
                                             parallel = TRUE,
                                             verbose = FALSE))

  pb$tick()



  # Yes, a for loop. The issue is that it adds repeatedly to the qseaGLM object, so can't be vectorised easily.
  print("Starting contrasts")
  for (i in seq_along(1:nrow(contrastsToDo))) {

    conName <- paste0(variable, contrastsToDo[i,"sample1"], "-", variable, contrastsToDo[i,"sample2"])

    if (!(contrastsToDo[i,"sample1"] %in% sampleTable[,variable])) {
      stop(glue::glue("value {contrastsToDo[i,]$sample1} not found in column {variable} of the sampleTable!"))
    }

    if (!(contrastsToDo[i,"sample2"] %in% sampleTable[,variable])) {
      stop(glue::glue("value {contrastsToDo[i,]$sample2} not found in column {variable} of the sampleTable!"))
    }

    if ("name" %in% colnames(contrastsToDo)) {
      conNameClean <- contrastsToDo[i,"name"] %>% dplyr::pull()
    }else{
      # Remove the hyphen from the name, because it messes up things later.
      # Also remove the variable name (whatever it is)
      conNameClean <- conName %>%
        stringr::str_replace("-","_vs_") %>%
        stringr::str_remove_all(variable)
    }

    print(glue::glue("Performing contrast {conNameClean}"))

    limContrast <- limma::makeContrasts(contrasts = conName, levels = design)

    qseaGLM <- suppressMessages(qsea::addContrast(qseaSet,
                                                  qseaGLM,
                                                  contrast = limContrast,
                                                  name = conNameClean,
                                                  parallel = TRUE,
                                                  verbose = FALSE))

    pb$tick()
    print("")
    if (mean(qseaGLM@contrast[[conNameClean]]$LRT_pval == 0) >= 0.1 & checkPVals) {

      if(is.null(covariates)){
        warning("Warning! More than 10% of windows have p-values of exactly 0, possibly something has gone wrong! \n
           Set checkPVals = FALSE to ignore this if sure.")
      } else {

      stop("Error! More than 10% of windows have p-values of exactly 0, possibly an error! \n
           Try not including covariates in the model if included, or set checkPVals = FALSE to ignore this if sure.")
      }
    }
  }
  print("Contrasts Complete")
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
#' @export
getDMRsData <- function(qseaSet, qseaGLM, sampleNames = NULL, variable = NULL, keepData = TRUE, keepGroupMeans = TRUE,
                        fdrThres = 0.05, keepPvals = FALSE, keepFragmentInfo = FALSE,
                        direction = "both"){

  #if (is.null(sampleTable)) {
  sampleTable <- qsea::getSampleTable(qseaSet)
  #}

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
  if (length(sigIndex) ==0 ) {
    print("No windows found, using hack to return, needs testing")
    sigIndex <- 1
    hack <- TRUE
  }

  if (length(sigIndex) ==1 ) {
    print("One windows found, using hack to return, needs testing")
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

  if(keepGroupMeans & variable != "group"){

    groupMeansList <- qsea::getSampleGroups(qseaSet)[names(qsea::getSampleGroups(qseaSet)) %in% unique(sampleTable$group)] %>%
      c(contrastMeansList)

  } else {
    groupMeansList <- contrastMeansList
  }


  dataTable <- qsea::makeTable(qseaSet, qseaGLM,
                               keep = sigIndex,
                               norm_methods = c("beta","nrpm"),
                               samples = sampleNames,
                               groupMeans = groupMeansList)

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

#' Make a full set of contrasts
#'
#' This function fits a GLM and returns the resulting data table
#'
#' @param qseaSet The qseaSet object
#' @param variable Which variable to use to calculate the DMRs between
#' @return A tibble with the data
#' @export
#'
makeAllContrasts <- function(qseaSet, variable){
  #TODO remove the message about new names by intercepting the
  contrastsToDo <- qseaSet %>%
    dropPooledControl() %>%
    qsea::getSampleTable() %>%
    tidyr::drop_na(variable) %>%
    dplyr::pull(variable) %>%
    unique() %>%
    gtools::mixedsort() %>%
    utils::combn(2) %>%
    t()  %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    stats::setNames(c("sample1","sample2")) %>%
    dplyr::filter(sample1 != sample2)

  return(contrastsToDo)
}

#' Fit the GLM and return the data
#'
#' This function fits a GLM and returns the resulting data table
#'
#' @param qseaSet The qseaSet object
#' @param variable Which variable to use to calculate the DMRs between
#' @param covariates Any variables to use as covariates, for instance patient in a paired analysis
#' @param contrastsToDo A data frame with two columns, sample1 and sample2, with the strings to compare between in each. Multiple rows means that multiple comparisons will be fitted
#' @param formula Alternative formula mode for calculating DMRs (not recommended)
#' @param minNRPM A minimum normalised reads per million value that at least one sample (in the contrasts) must reach in order to consider the region for further calculation. Set this or minReadCount (but preferably this).
#' @param minReadCount A minimum read count that at least one sample (in the contrasts) must reach in order to consider the region for further calculation. Preferably use minNRPM.
#' @param sampleNames Which samples to return the information for.
#' @param fdrThres False discovery rate threshold.
#' @param keepData Whether to keep the individual data columns in the output
#' @param keepGroupMeans Whether to keep the group means in the output
#' @param keepPvals Whether to keep the unadjusted p-values in the output
#' @param checkPVals Whether to check that the p-values aren't all zero to avoid a bug with covariates, only turn this off if you are very sure what you are doing!
#' @param direction Whether to keep regions that are up/down/both.
#' @return A tibble with the data
#' @export
#'
calculateDMRs <- function(qseaSet, variable = NULL,
                          sampleNames = NULL,
                          covariates = NULL,
                          contrastsToDo = NULL,
                          minReadCount = 0,
                          minNRPM = 1,
                          checkPVals = TRUE,
                          fdrThres = 0.05,
                          keepPvals = FALSE,
                          formula = NULL,
                          keepData = FALSE,
                          keepGroupMeans = FALSE,
                          direction = "both"){

  if (is.null(variable)) {stop("variable must be specified!")}

  if (!is.data.frame(contrastsToDo)) {
    if (contrastsToDo %in% c("All", "all")) {
      contrastsToDo <- makeAllContrasts(qseaSet, variable)
      message(
        glue::glue(
          "Calculating {nrow(contrastsToDo)} possible contrasts on the {variable} column."
        )
      )
    } else if (contrastsToDo %in% c("First", "first")) {
        contrastsToDo <- makeAllContrasts(qseaSet, variable)[1, ]
        message(glue::glue(
          "Calculating the first possible contrast on the {variable} column."
        ))
    } else if (stringr::str_detect(contrastsToDo,"All_vs_")){

        value2 = contrastsToDo %>% stringr::str_remove("All_vs_")
        contrastsToDo <- tibble::tibble(sample1 = qseaSet %>% pullQset(variable) %>% unique() %>% setdiff(value2),
                                        sample2 = value2)
    }

  }



  if(is.null(contrastsToDo)){stop("No contrasts specified!")}

  qseaGLM <- fitQseaGLM(qseaSet, variable = variable,  covariates = covariates,
                        contrastsToDo = contrastsToDo,  minReadCount = minReadCount,
                        minNRPM = minNRPM,
                        checkPVals = checkPVals, formula = formula)

  dataTable <- getDMRsData(qseaSet, qseaGLM, sampleNames = sampleNames,
                           fdrThres = fdrThres, keepPvals = keepPvals, keepData = keepData,
                           keepGroupMeans = keepGroupMeans,
                           variable = variable,
                           direction = direction) %>%
    dplyr::rename(seqnames = chr, start = window_start, end = window_end)

  deltas  <- purrr::map_dfc(1:nrow(contrastsToDo),
                            function(x){
                              name1 <- contrastsToDo[x,]$sample1
                              name2 <- contrastsToDo[x,]$sample2
                              (dataTable[,paste0(name1,"_beta_means")] - dataTable[,paste0(name2,"_beta_means")]) %>%
                                tibble::enframe(name = "rowIndex", value = paste0(name1,"_vs_",name2,"_betaDelta")) %>%
                                dplyr::select(-rowIndex)
                            }
  )

  dataTable %>%
    dplyr::bind_cols(deltas) %>%
    tibble::as_tibble() %>%
    return()

}
