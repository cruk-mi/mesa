#' Fit a generalized linear model (GLM) to a qseaSet
#'
#' Fit a negative-binomial GLM using [qsea::fitNBglm()] on counts from a
#' [qsea::qseaSet], and add one or more contrasts.  
#' This is the low-level worker used by [calculateDMRs()].
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing counts and sample metadata.
#'
#' @param variable `character(1)`  
#'   Column name in the sample table used as the primary explanatory variable.  
#'   **Default:** `NULL` (must be provided unless `formula` is supplied).
#'
#' @param covariates `character()`  
#'   Additional column names (e.g., batch, patient ID) to include as covariates.  
#'   **Default:** `NULL`.
#'
#' @param formula `formula` or `NULL`  
#'   Optional model formula overriding `variable`/`covariates`.  
#'   **Default:** `NULL`.
#'
#' @param contrasts `data.frame` or `NULL`  
#'   A data frame with columns `group1` and `group2` specifying contrasts to
#'   compute. If `NULL`, a fitted GLM with no contrasts added is returned.  
#'   **Default:** `NULL`.
#'
#' @param keepIndex `integer()` or `NULL`  
#'   Row indices (windows) to retain **before** fitting; when supplied, this
#'   overrides filtering by `minReadCount` / `minNRPM`.  
#'   **Default:** `NULL`.
#'
#' @param minReadCount `numeric(1)`  
#'   Minimum raw count required in at least one sample to retain a window.  
#'   **Default:** `0`.
#'
#' @param minNRPM `numeric(1)`  
#'   Minimum NRPM required in at least one sample to retain a window.  
#'   **Default:** `1`.
#'
#' @param checkPVals `logical(1)`  
#'   If `TRUE`, stop or warn if an excessive proportion of p-values are exactly
#'   zero (possible instability).  
#'   **Default:** `TRUE`.
#'
#' @param calcDispersionAll `logical(1)`  
#'   If `TRUE`, samples not present in any specified contrast are still used to
#'   fit the initial GLM and dispersion estimates. Note this means adding samples
#'   to `qseaSet` can change DMRs even when they are not directly contrasted.  
#'   **Default:** `FALSE`.
#'
#' @return A `qseaGLM` object containing:
#' * the fitted negative-binomial GLM;  
#' * (optionally) one or more added contrasts (`group1` vs `group2`);  
#' * filtering metadata (based on `keepIndex`, `minReadCount`, `minNRPM`);  
#' * model specification (`variable`, `covariates`, or `formula`).  
#'
#' @details
#' * Dispersion is estimated across the windows retained after filtering.  
#' * When `calcDispersionAll = TRUE`, all samples contribute to dispersion
#'   estimation even if not part of any contrast.  
#' * Contrasts are added sequentially with [qsea::addContrast()].  
#' * By default, beta-normalised counts are used as the outcome (as in qsea).
#'
#' @family DMR-detection
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' contr <- tibble::tibble(group1 = "LUAD", group2 = "NormalLung")
#' fitQseaGLM(qs, variable = "type", contrasts = contr)
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
        warning("More than 20% of windows have p-values of exactly 0, possibly something has gone wrong! \n
           Set checkPVals = FALSE to ignore this.")
      } else {

        stop(
          "More than 20% of windows have p-values of exactly 0; this likely ",
          "indicates a model issue.\n",
          "Try removing covariates from the model (if any), or set ",
          "checkPVals = FALSE to ignore this if you're sure."
        )
      }
    }
  }
  message("Contrasts Complete")
  return(qseaGLM)
}


#' Extract DMR-level results from a fitted GLM
#'
#' Given a [`qseaGLM`] object with one or more contrasts, extract a *wide*
#' results table of statistics, adjusted p-values, and (optionally) group means
#' and per-sample values.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object used to obtain metadata and counts.
#'
#' @param qseaGLM `qseaGLM`.  
#'   Fitted model returned by [fitQseaGLM()].
#'
#' @param sampleNames `character()` or `NULL`.  
#'   Sample names to include as per-sample columns.  
#'   **Default:** `NULL` (chosen automatically when `keepData = TRUE`, otherwise
#'   per-sample columns are omitted).
#'
#' @param variable `character(1)` or `NULL`.  
#'   Sample-table column on which contrasts were based (used to compute group
#'   means).  
#'   **Default:** `NULL` (required if `keepGroupMeans = TRUE`).
#'
#' @param FDRthres `numeric(1)`.  
#'   False discovery rate threshold used to flag significance.  
#'   **Default:** `0.05`.
#'
#' @param keepData `logical(1)`.  
#'   If `TRUE`, include per-sample values (beta/NRPM etc.) as columns.  
#'   **Default:** `FALSE`.
#'
#' @param keepGroupMeans `logical(1)`.  
#'   If `TRUE`, include group mean columns for each contrast side (based on
#'   `variable`).  
#'   **Default:** `FALSE`.
#'
#' @param keepPvals `logical(1)`.  
#'   If `TRUE`, include raw (unadjusted) p-values in addition to FDR.  
#'   **Default:** `FALSE`.
#'
#' @param keepFragmentInfo `logical(1)`.  
#'   If `TRUE`, include fragment/MAPQ metrics where available.  
#'   **Default:** `FALSE`.
#'
#' @param direction `character(1)`.  
#'   Direction for significance calling, passed to [qsea::isSignificant()]:
#'   one of `"up"`, `"down"`, or `"both"`.  
#'   **Default:** `"both"`.
#'
#' @return A tibble with one row per DMR (window) containing:
#' * **Per-contrast statistics** (e.g., log2FC, deltaBeta, test statistic).  
#' * **Multiple-testing output** (FDR; optionally raw p-values if `keepPvals`).  
#' * **Optional group summaries** if `keepGroupMeans = TRUE`.  
#' * **Optional per-sample columns** if `keepData = TRUE` (for `sampleNames`).  
#' * **Optional fragment/MAPQ metrics** if `keepFragmentInfo = TRUE`.  
#'
#' @details
#' Significant rows are determined using `FDRthres` and `direction` via
#' [qsea::isSignificant()]. When `keepData = TRUE` and `sampleNames = NULL`,
#' per-sample columns are included for all samples in the fitted model (or the
#' subset relevant to each contrast, depending on the internal representation).
#' Group means require `variable` to identify the groupings used in contrasts.
#'
#' @family DMR-detection
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' contr <- tibble::tibble(group1 = "LUAD", group2 = "NormalLung")
#' glmfit <- fitQseaGLM(qs, variable = "type", contrasts = contr)
#' getDMRsData(qs, glmfit, variable = "type", FDRthres = 0.1)
#' 
#' @export
getDMRsData <- function(qseaSet, qseaGLM, sampleNames = NULL, variable = NULL, keepData = FALSE, keepGroupMeans = FALSE,
                        FDRthres = 0.05, keepPvals = FALSE, keepFragmentInfo = FALSE,
                        direction = "both"){

  sampleTable <- qsea::getSampleTable(qseaSet)

  if (is.null(sampleNames)){
    if(keepData) {
      sampleNames <- sampleTable$sample_name
    } else {
      sampleNames <- NULL
    }
  }

  sigIndex <- purrr::map(names(qseaGLM@contrast), ~ qsea::isSignificant(qseaGLM, contrast = ., fdr_th = FDRthres, direction = direction)) %>%
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
#' @param qseaSet `qseaSet`.  
#'   Input object providing the sample table from which factor levels are taken.
#'
#' @param variable `character(1)`.  
#'   Column name in the sample table whose distinct levels define the contrasts.  
#'   **Default:** none (must be supplied).
#'
#' @return A tibble with two columns:  
#' * `group1` – the first level of the pair,  
#' * `group2` – the second level of the pair.  
#' Each row represents one **unordered** pair of distinct levels (i.e., all
#' combinations without repetition).
#'
#' @details
#' Levels are inferred from `unique(qsea::getSampleTable(qseaSet)[[variable]])`.
#' Pairs are formed between all distinct levels; no self-contrasts are included.
#'
#' @family DMR-detection
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # All pairwise contrasts among levels of 'type'
#' exampleTumourNormal %>%
#'   makeAllContrasts("type")
#'
#' @export
makeAllContrasts <- function(qseaSet, variable){
  vals <- qseaSet %>%
    qsea::getSampleTable() %>%
    filter(!is.na({{variable}})) %>%
    dplyr::pull({{variable}}) %>%
    unique() %>%
    gtools::mixedsort()

  tidyr::expand_grid(group1 = vals, group2 = vals) %>%
    filter(group1 < group2) %>%
    return()
}


#' Fit GLM and return DMR data in one step
#'
#' High-level wrapper that calls [fitQseaGLM()] and [getDMRsData()] to fit a
#' negative-binomial GLM and extract a DMR results table in one go.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing counts and sample metadata.
#'
#' @param variable `character(1)` or `NULL`.  
#'   Sample-table column used as the primary explanatory variable for contrasts.  
#'   **Default:** `NULL` (must be supplied unless `formula` is provided).
#'
#' @param covariates `character()` or `NULL`.  
#'   Additional covariate column names (e.g., batch, patient).  
#'   **Default:** `NULL`.
#'
#' @param contrasts `data.frame`, `character(1)`, or `NULL`.  
#'   Contrast specification. If a data frame, must contain columns `group1` and
#'   `group2`. If a string, supports:
#'   * `"A_vs_B"` — a single explicit contrast,
#'   * `"All"` / `"all"` — all pairwise contrasts among levels of `variable`,
#'   * `"All_vs_X"` — each level vs `X`,
#'   * `"X_vs_All"` — `X` vs each other level,
#'   * `"first"` — only the first possible contrast (based on factor order).  
#'   **Default:** `NULL` (fit the GLM without adding contrasts).
#'
#' @param minReadCount `numeric(1)`.  
#'   Minimum raw count required in at least one sample to keep a window.  
#'   **Default:** `0`.
#'
#' @param minNRPM `numeric(1)`.  
#'   Minimum NRPM required in at least one sample to keep a window.  
#'   **Default:** `1`.
#'
#' @param checkPVals `logical(1)`.  
#'   If `TRUE`, stop/warn when an excessive proportion of p-values are exactly
#'   zero (instability check).  
#'   **Default:** `TRUE`.
#'
#' @param FDRthres `numeric(1)`.  
#'   False discovery rate threshold for calling significance.  
#'   **Default:** `0.05`.
#'
#' @param keepPvals `logical(1)`.  
#'   If `TRUE`, include raw (unadjusted) p-values in the output.  
#'   **Default:** `FALSE`.
#'
#' @param formula `formula` or `NULL`.  
#'   Optional model formula overriding `variable`/`covariates`.  
#'   **Default:** `NULL`.
#'
#' @param keepContrastMeans `logical(1)`.  
#'   If `FALSE`, drop contrast mean columns from the output.  
#'   **Default:** `TRUE`.
#'
#' @param keepData `logical(1)`.  
#'   If `TRUE`, include per-sample data columns for significant windows.  
#'   **Default:** `FALSE`.
#'
#' @param keepGroupMeans `logical(1)`.  
#'   If `TRUE`, include group mean columns (based on `variable`).  
#'   **Default:** `FALSE`.
#'
#' @param direction `character(1)`.  
#'   Direction for calling significance, passed to [qsea::isSignificant()]:
#'   one of `"up"`, `"down"`, or `"both"`.  
#'   **Default:** `"both"`.
#'
#' @param calcDispersionAll `logical(1)`.  
#'   If `TRUE`, samples not present in any contrast still contribute to the
#'   initial GLM fit/dispersion estimation (so adding samples can change DMRs).  
#'   **Default:** `FALSE`.
#'
#' @return A tibble with one row per **significant** region/window containing:
#' * per-contrast statistics (e.g., log2FC, deltaBeta, test statistic),  
#' * multiple-testing results (FDR; optionally raw p-values if `keepPvals`),  
#' * optional group means (if `keepGroupMeans = TRUE`),  
#' * optional per-sample columns (if `keepData = TRUE`),  
#' * metadata columns for genomic coordinates and window attributes.
#'
#' @details
#' Contrast strings are parsed as follows:
#' * `"A_vs_B"` creates a single contrast *A vs B*;  
#' * `"All"`/`"all"` expands to all pairwise contrasts among levels of `variable`
#'   (via [makeAllContrasts()]);  
#' * `"All_vs_X"` creates *each level vs X*;  
#' * `"X_vs_All"` creates *X vs each other level*;  
#' * `"first"` constructs only the first possible pair (based on factor order).  
#' Internally, this function:
#' 1) filters windows by `minReadCount`/`minNRPM` (unless overridden inside the
#'    GLM wrapper), 2) fits a NB-GLM with [fitQseaGLM()], 3) extracts a wide
#'    DMR table with [getDMRsData()], and 4) prunes contrast means if
#'    `keepContrastMeans = FALSE`.
#'
#' @family DMR-detection
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # A single explicit contrast
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type", contrasts = "LUAD_vs_NormalLung")
#'
#' # All pairwise contrasts among levels of 'type'
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type", contrasts = "all")
#'
#' # Tighter FDR threshold
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type", contrasts = "all", FDRthres = 1e-4)
#'
#' # Return all per-sample columns for the first available contrast on 'tumour'
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "tumour", contrasts = "first", keepData = TRUE)
#'
#' @export
calculateDMRs <- function(qseaSet,
                          variable = NULL,
                          covariates = NULL,
                          contrasts = NULL,
                          minReadCount = 0,
                          minNRPM = 1,
                          checkPVals = TRUE,
                          FDRthres = 0.05,
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
                           FDRthres = FDRthres, keepPvals = keepPvals, keepData = keepData,
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
