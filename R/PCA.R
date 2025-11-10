#' Generate a PCA from a qseaSet
#'
#' Convenience wrapper around [getDimRed()] with `method = "PCA"`.
#'
#' @describeIn getDimRed Principal component analysis of a qseaSet.
#' @family dimred-helpers
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Default PCA (beta-normalised, top 1000 most variable windows, 5 PCs)
#' exampleTumourNormal %>%
#'   getPCA()
#'
#' # Request only the top 3 PCs
#' exampleTumourNormal %>%
#'   getPCA(nPC = 3)
#'
#' # Use multiple cutoffs for most-variable windows (returns list of PCA results)
#' exampleTumourNormal %>%
#'   getPCA(topVarNum = c(10, 100, 500, NA)) %>%
#'   slot("res") %>%
#'   names()
#' @export
getPCA <- function(qseaSet,
                   dataTable = NULL,
                   regionsToOverlap = NULL,
                   normMethod = "beta",
                   minEnrichment = 3,
                   useGroupMeans = FALSE,
                   minDensity = 0,
                   topVarNum = 1000,
                   topVarSamples = NULL,
                   center = TRUE,
                   scale = FALSE,
                   nPC = 5,
                   returnDataTable = FALSE, 
                   verbose = TRUE) {

  getDimRed(qseaSet = qseaSet,
            dataTable = dataTable,
            method = "PCA",
            regionsToOverlap = regionsToOverlap,
            normMethod = normMethod,
            minEnrichment = minEnrichment,
            useGroupMeans = useGroupMeans,
            minDensity = minDensity,
            topVarNum = topVarNum,
            topVarSamples = topVarSamples,
            center = center,
            scale = scale,
            nPC = nPC,
            returnDataTable = returnDataTable,
            verbose = verbose)

}


#' Generate a UMAP from a qseaSet
#'
#' Convenience wrapper around [getDimRed()] with `method = "UMAP"`.
#'
#' @describeIn getDimRed Uniform manifold approximation and projection of a qseaSet.
#' @family dimred-helpers
#' @examples
#' # Quick demo on a small synthetic qseaSet
#' qsea::getExampleQseaSet(repl = 20) %>%
#'   getUMAP()
#'
#' # Alter UMAP hyperparameters
#' qsea::getExampleQseaSet(repl = 20) %>%
#'   getUMAP(n_neighbors = 5, min_dist = 1)
#'
#' # Use top 500 most variable windows only
#' qsea::getExampleQseaSet(repl = 20) %>%
#'   getUMAP(topVarNum = 500)
#' @export
getUMAP <- function(qseaSet,
                   dataTable = NULL,
                   regionsToOverlap = NULL,
                   normMethod = "beta",
                   minEnrichment = 3,
                   useGroupMeans = FALSE,
                   minDensity = 0,
                   topVarNum = 1000,
                   topVarSamples = NULL,
                   returnDataTable = FALSE,
                   verbose = TRUE,
                   ...) {

  getDimRed(qseaSet = qseaSet,
            dataTable = dataTable,
            method = "UMAP",
            regionsToOverlap = regionsToOverlap,
            normMethod = normMethod,
            minEnrichment = minEnrichment,
            useGroupMeans = useGroupMeans,
            minDensity = minDensity,
            topVarNum = topVarNum,
            topVarSamples = topVarSamples,
            returnDataTable = returnDataTable,
            ...)

}


#' Dimensionality reduction for qseaSets
#'
#' Run PCA or UMAP on methylation signal matrices extracted from a
#' [qsea::qseaSet]. This is the core workhorse for dimensionality reduction;
#' wrappers such as [getPCA()] and [getUMAP()] provide simplified access.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object providing windows, counts and `sampleTable`.
#'
#' @param dataTable `data.frame` or `GRanges` or `NULL`.  
#'   Normalised values with windows in rows and samples in columns. Must contain
#'   `seqnames`, `start`, `end` (if `data.frame`), or metadata columns with
#'   values (if `GRanges`). If `NULL`, the matrix is derived internally (see
#'   `normMethod`, `minEnrichment`).  
#'   **Default:** `NULL`.
#'
#' @param method `character(1)`. One of `"PCA"` or `"UMAP"`.  
#'   **Default:** `"PCA"`.
#'
#' @param regionsToOverlap `GRanges`, coercible `data.frame`, or `NULL`.  
#'   If supplied, only windows overlapping these regions are used.  
#'   **Default:** `NULL`.
#'
#' @param normMethod `character(1)`.  
#'   Name of predefined normalisation (e.g., `"beta"`, `"nrpm"`); passed to
#'   [qsea::normMethod()] when `dataTable = NULL`.  
#'   **Default:** `"beta"`.
#'
#' @param minEnrichment `numeric(1)`.  
#'   Minimum reads required for beta values to be non-`NA`; forwarded to
#'   [getDataTable()] when building data internally.  
#'   **Default:** `3`.
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, average samples within the `group` column (combine replicates)
#'   before DR.  
#'   **Default:** `FALSE`.
#'
#' @param minDensity `numeric(1)`.  
#'   Minimum CpG density; windows below this are removed.  
#'   **Default:** `0`.
#'
#' @param topVarNum `numeric(1)` or `integer(1)` or `Inf` or `NULL` **or** a vector.  
#'   Keep the `topVarNum` most variable windows (by SD). If `NA`, `NULL`, `Inf`,
#'   or `>=` available windows, use all. If a vector, DR is run once per value
#'   and results are returned in a list.  
#'   **Default:** `1000`.
#'
#' @param topVarSamples `NULL`, `character()`, `list`, or `regex` string.  
#'   Samples used to compute variability. If `NULL`/`NA`, use all samples. If a
#'   vector/regex, select matching samples. If a `list`, perform DR per list
#'   element (should match `length(topVarNum)` when vectorised).  
#'   **Default:** `NULL`.
#'
#' @param center `logical(1)`.  
#'   Centre features to mean zero (PCA only).  
#'   **Default:** `TRUE`.
#'
#' @param scale `logical(1)`.  
#'   Scale features to unit variance (PCA only).  
#'   **Default:** `FALSE`.
#'
#' @param nPC `integer(1)`.  
#'   Number of principal components to compute (PCA only).  
#'   **Default:** `5`.
#'
#' @param returnDataTable `logical(1)`.  
#'   If `TRUE`, include the matrix used for DR in the return object.  
#'   **Default:** `FALSE`.
#'   
#' @param verbose `logical(1)`.  
#'   If `TRUE`, print messages regarding the function execution.  
#'   **Default:** `TRUE`.
#'   
#' @param ... Additional arguments passed to [uwot::umap()] (when `method="UMAP"`),
#'   e.g. `n_neighbors`, `min_dist`, `metric`.
#'
#' @return A [`mesaDimRed`] object with:
#' * **res**: Named list of DR results; each element is a [`mesaPCA`] (for PCA)
#'   or [`mesaUMAP`] (for UMAP). Names correspond to entries in `params$topVar`.
#' * **samples**: Character vector of sample IDs used.
#' * **dataTable** (optional): The numeric matrix used for DR when
#'   `returnDataTable = TRUE`. It reflects filtering by `regionsToOverlap` and
#'   `minDensity`, removal of rows with missing values, and may include columns
#'   of window SDs when applicable (named per `params$topVar`).
#' * **params**: List of parameters used (method, selection/filter settings, etc.).
#' * **sampleTable**: Copy of `qsea::getSampleTable(qseaSet)`.
#'
#' @details:
#' - For `method = "PCA"`, DR is performed via [stats::prcomp()] with `center`
#'   and `scale` controls.  
#' - For `method = "UMAP"`, DR is performed via [uwot::umap()].  
#' - Prior to DR, optional filtering is applied: enrichment cut-off
#'   (`minEnrichment`), CpG density (`minDensity`), region restriction
#'   (`regionsToOverlap`), and most-variable windows (`topVarNum`).  
#' - When `topVarNum` or `topVarSamples` are vectorised/lists, multiple DR runs
#'   are performed and collected in `res`.
#'
#' @seealso
#' [getPCA()], [getUMAP()], [mesaDimRed-class], [mesaPCA-class], [mesaUMAP-class],
#' [qsea::normMethod()], [uwot::umap()], [stats::prcomp()]
#'
#' @family dimred-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # PCA on default beta-normalised matrix, keep 3 PCs
#' exampleTumourNormal %>%
#'   getDimRed(method = "PCA", nPC = 3)
#'
#' # UMAP on the top 500 most variable windows
#' exampleTumourNormal %>%
#'   getDimRed(method = "UMAP", topVarNum = 500, n_neighbors = 5)
#'
#' # Restrict DR to a set of regions (coerced from data.frame)
#' regs <- data.frame(seqnames = "chr7", start = 25002001, end = 25027900)
#' exampleTumourNormal %>%
#'   getDimRed(method = "PCA", regionsToOverlap = regs)
#'
#' # Vectorised topVarNum runs (returns a list of results)
#' exampleTumourNormal %>%
#'   getDimRed(method = "UMAP", topVarNum = c(500, 2000), n_neighbors = 5) %>%
#'   slot("res") %>%
#'   names()
#' @export
getDimRed <- function(qseaSet,
                      dataTable = NULL,
                      method = "PCA",
                      regionsToOverlap = NULL,
                      normMethod = "beta",
                      minEnrichment = 3,
                      useGroupMeans = FALSE,
                      minDensity = 0,
                      topVarNum = 1000,
                      topVarSamples = NULL,
                      center = TRUE,
                      scale = FALSE,
                      nPC = 5,
                      returnDataTable = FALSE,
                      verbose = TRUE,
                      ...) {

  if (!("qseaSet" %in% class(qseaSet))) {
    stop("Please provide a qseaSet in the first position.")
  }

  if (qseaSet %>% qsea::getSampleNames() %>% length() <= 2) {
    stop(glue::glue("Insufficient samples {qseaSet %>% qsea::getSampleNames() %>% length()} in provided qseaSet, must be at least 3."))
  }

  if (useGroupMeans) {
    groupString <- "sample group"
  } else {
    groupString <- "sample"
  }

  if (!is.null(dataTable)) {
    dataTable <- asValidGranges(dataTable)

    samples <- dataTable %>%
      tidyr::as_tibble() %>%
      dplyr::select(-tidyselect::any_of(c("seqnames", "start", "end", "width", "strand", "CpG_density"))) %>%
      colnames()

    normMethodSuffixDetected <- stringr::str_detect(samples, glue::glue("_{normMethod}$"))

    if (all(normMethodSuffixDetected)) {
      samples <- stringr::str_remove(samples, glue::glue("_{normMethod}$"))
      dataTable <- removeNormMethodSuffix(dataTable %>% tidyr::as_tibble(),
                                          normMethod) %>%
        asValidGranges()
    } else if (any(normMethodSuffixDetected)) {
      stop(glue::glue("normMethod suffix '_{normMethod}' is present for some but not all of the {groupString} column names in dataTable."))
    }

    if (useGroupMeans) {
      if (!all(samples %in% names(qsea::getSampleGroups(qseaSet)))) {
        stop(glue::glue("At least one {groupString} name in dataTable does not have a matching {groupString} name in qseaSet.
             (If there is a normMethod suffix on the {groupString} column names in dataTable, check it matches the input normMethod argument: '_{normMethod}')"))
      }
    } else {
      if (!all(samples %in% qsea::getSampleNames(qseaSet))) {
        stop(glue::glue("At least one {groupString} name in dataTable does not have a matching {groupString} name in qseaSet.
             (If there is a normMethod suffix on the {groupString} column names in dataTable, check it matches the input normMethod argument: '_{normMethod}')"))
      }
    }

    if (length(plyranges::setdiff_ranges(dataTable, getWindows(qseaSet))) > 0) {
      stop("At least one window in dataTable does not have a matching window in qseaSet.")
    }

    testSamples <- samples[1:min(5, length(samples))]
    inputValues <- dataTable[1:min(50, length(dataTable)), ]
    testValues <- getDataTable(qseaSet %>%
                                 filter(sample_name %in% testSamples) %>%
                                 filterByOverlaps(inputValues),
                               normMethod,
                               useGroupMeans = useGroupMeans) %>%
      dplyr::arrange(seqnames, start, end)

    if (!isTRUE(all.equal(testValues, inputValues %>%
                          tidyr::as_tibble() %>%
                          dplyr::select(tidyselect::all_of(colnames(testValues))) %>%
                          dplyr::arrange(seqnames, start, end)))) {
      warning("Newly-generated dataTable from qseaSet on a small subset of samples/windows does not match input dataTable.")
    }
    initialNumWindows <- length(dataTable)
  } else {

    if (useGroupMeans) {
      samples <- names(qsea::getSampleGroups(qseaSet))
    } else {
      samples <- qsea::getSampleNames(qseaSet)
    }
    initialNumWindows <- getWindows(qseaSet) %>% length()
  }

  if(verbose){
  message(glue::glue("------------------------------
                     Initial number of windows = {initialNumWindows}."))
  }

  if (!is.null(regionsToOverlap)) {
    regionsToOverlap <- regionsToOverlap %>%
      tibble::as_tibble() %>%
      dplyr::select(tidyselect::any_of(c("seqnames", "start", "end", "CpG_density"))) %>%  # keep only minimum columns necessary
      asValidGranges()

    if (is.null(dataTable)) {
      qseaSet <- qseaSet %>%
        filterByOverlaps(regionsToOverlap = regionsToOverlap)

      numWindowsRemovedRegionOverlap <- initialNumWindows - length(getWindows(qseaSet))
      if(verbose) {
        message(glue::glue("Filtered out {numWindowsRemovedRegionOverlap} windows using regionsToOverlap: {length(getRegions(qseaSet))} windows remaining."))
      }

        
    } else {
      dataTable <- dataTable %>%
        plyranges::filter_by_overlaps(y = regionsToOverlap)

      numWindowsRemovedRegionOverlap <- initialNumWindows - length(dataTable)
      if(verbose) {
        message(glue::glue("Filtered out {numWindowsRemovedRegionOverlap} windows using regionsToOverlap: {length(dataTable)} windows remaining."))
      }
    }
  } else {
    numWindowsRemovedRegionOverlap <- NULL
  }

  if (is.null(dataTable)) {
    if(verbose) { message("-----------") }
    dataTable <- qseaSet %>%
      filterWindows(CpG_density >= minDensity) %>%
      getDataTable(normMethod = normMethod, useGroupMeans = useGroupMeans)
    if(verbose) { message("-----------") }

    if (minDensity > 0) {
      numWindowsRemovedMinDensity <- length(getWindows(qseaSet)) - nrow(dataTable)
      if(verbose) {
        message(glue::glue("Filtered out {numWindowsRemovedMinDensity} windows with CpG_density < {minDensity}: {nrow(dataTable)} windows remaining."))
      }
    } else {
      numWindowsRemovedMinDensity <- NULL
    }

  } else {
    currentNumWindows <- length(dataTable)
    dataTable <- dataTable %>%
      tidyr::as_tibble() %>%
      dplyr::left_join(qseaSet %>%
                         getWindows() %>%
                         tidyr::as_tibble() %>%
                         dplyr::select(seqnames, start, end, CpG_density),
                       by = dplyr::join_by(seqnames, start, end, CpG_density)) %>%
      dplyr::filter(CpG_density >= minDensity)

    if (minDensity > 0) {
      numWindowsRemovedMinDensity <- currentNumWindows - nrow(dataTable)
      if(verbose) {
        message(glue::glue("Filtered out {numWindowsRemovedMinDensity} windows with CpG_density < {minDensity}: {nrow(dataTable)} windows remaining."))
      }
    } else {
      numWindowsRemovedMinDensity <- NULL
    }

  }

  currentNumWindows <- nrow(dataTable)

  dataTable <- dataTable %>%
    tidyr::drop_na(tidyr::all_of(samples))

  numWindowsRemovedMissingVals <- currentNumWindows - nrow(dataTable)

  if (numWindowsRemovedMissingVals > 0) {
    if(verbose){
      message(glue::glue("Filtered out {numWindowsRemovedMissingVals} windows with at least one missing value: {nrow(dataTable)} windows remaining.
                       ------------------------------"))
    }
  } else {
    if(verbose){
      message(glue::glue("No windows have missing values.
            ------------------------------"))
    }
  }

  if (nrow(dataTable) <= 2) {
    stop(glue::glue("Insufficient windows {nrow(dataTable)} remaining after filtering! Have you filtered for poor quality samples?"))
  }

  if (is.null(topVarNum)) {
    topVarNum <- NA
  } else if (!(length(topVarNum) == 1 && is.na(topVarNum)) & !is.vector(topVarNum, mode = "numeric")) {
    stop("topVarNum should be a numeric vector (or NULL")
  }

  if (!is.list(topVarSamples)) {
    topVarSamples <- list(topVarSamples)
  }

  if (length(topVarNum) > 1 && !(length(topVarSamples) %in% c(1, length(topVarNum)))) {
    stop("If topVarSamples is a list and length(topVarNum) > 1, topVarSamples should be the same length as topVarNum.")
  }


  # replace NA with NULL
  topVarSamples <- purrr::map(topVarSamples, function(tVS) {

    if (length(tVS) == 1 && is.na(tVS)) {
      tVS <- NULL
    }

    return(tVS)
  })

  topVarSamplesInput <- topVarSamples

  topVarSamples <- purrr::map(topVarSamples, function(tVS) {

    if (!is.null(tVS)) {
      if (!is.vector(tVS, mode = "character")) {
        stop(glue::glue("topVarSamples should be NULL, a character vector of {groupString} names or a regular expression (or a list of these)."))

      } else if (length(tVS) > 1) { # character vector of sample (group) names
        notInSamples <- setdiff(tVS, samples)
        if (length(notInSamples) > 0) {
          stop(glue::glue("topVarSamples contains {groupString} names that are not in the qseaSet and/or dataTable:
                        {paste0(notInSamples, collapse = ', ')}."))
        }

      } else { # regular expression to match
        tVS <- stringr::str_subset(samples, tVS)

      }

    } else {
      tVS <- samples

    }

    return(tVS)

  })

  topVar <- tibble::tibble(topVarNum = topVarNum, topVarSamples = topVarSamples, topVarNumInput = topVarNum, topVarSamplesInput = topVarSamplesInput) %>%
    dplyr::mutate(topVarNum = ifelse(topVarNum >= nrow(dataTable), NA, topVarNum),
                  topVarSamples = ifelse(is.na(topVarNum), list(NULL), topVarSamples),
                  inputChanged = !purrr::map2_lgl(topVarNum, topVarNumInput, identical) |
                    !purrr::map2_lgl(topVarSamples, topVarSamplesInput, identical))

  if (any(topVar$topVarNumInput > nrow(dataTable), na.rm = TRUE)) {
    message(glue::glue("The following topVarNum values are larger than the number of remaining windows (= {nrow(dataTable)}): {paste0(dplyr::filter(topVar, topVarNumInput > nrow(dataTable)) %>% pull(topVarNumInput) %>% unique(), collapse = ', ')}"))

    if (any(is.na(topVar$topVarNumInput) | topVar$topVarNumInput == nrow(dataTable))) {
      message(glue::glue("These values are not used; {method} is already being done with all remaining windows."))
    } else {
      message(glue::glue("These values are not used; {method} will be done with all remaining windows instead."))
    }
  }

  topVar <- topVar %>%
    dplyr::distinct(topVarNum, topVarSamples, .keep_all = TRUE) %>%
    dplyr::group_by(topVarSamples) %>%
    dplyr::mutate(windowSdName = ifelse(!is.na(topVarNum), glue::glue("windowSd{dplyr::cur_group_id()}"), NA), .before = 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(windowSdName, topVarNum) %>%
    dplyr::mutate(resName = glue::glue("{method %>% stringr::str_to_lower()}{dplyr::row_number()}"), .before = 1) 

  rowSds <- function(x) {
    sqrt(rowSums((x - rowMeans(x)) ^ 2) / (ncol(x) - 1))
  }

  if (any(!is.na(topVar$windowSdName))) {
    dataTable <- topVar %>%
      tidyr::drop_na(windowSdName) %>%
      dplyr::select(topVarSamples, windowSdName) %>%
      dplyr::distinct() %>%
      dplyr::mutate(topVarSamples = purrr::set_names(topVarSamples, windowSdName)) %>%
      dplyr::pull(topVarSamples) %>%
      purrr::imap(function(tVS, nm) {

        if (length(setdiff(samples, tVS)) == 0) {
          if (verbose) {
            message(glue::glue("Calculating standard deviation for each window across all {length(samples)} {groupString}s:
                           {paste0(tVS, collapse = ', ')}.
                           -> column name {nm}."))
          }
        } else {
          if(verbose) {
            message(glue::glue("Calculating standard deviation for each window across {length(tVS)} of {length(samples)} {groupString}s:
                           {paste0(tVS, collapse = ', ')}.
                           -> column name {nm}."))
          }
        }

        if (length(tVS) <= 2) {
          stop(glue::glue("Standard deviation calculated on less than 3 {groupString} names. Insufficent number of {groupString}s to calculate variance."))
        }

        dataTable <- dataTable %>%
          dplyr::mutate({{nm}} := dplyr::select(., tidyr::all_of(tVS)) %>%
                          rowSds()) %>%
          dplyr::select(seqnames, start, end, tidyr::all_of(nm))

      }) %>%
      purrr::reduce(dplyr::left_join) %>%
      dplyr::left_join(dataTable, .,
                       by = dplyr::join_by(seqnames, start, end))
  }

  res <- topVar %>%
    dplyr::group_by(windowSdName) %>%
    dplyr::group_map(.keep = TRUE, .f = function(sdGp, sdName) {

      if (!is.na(sdName$windowSdName)) {
        dataTable <- dataTable %>%
          dplyr::arrange(dplyr::desc(!!dplyr::sym(sdName$windowSdName)))
      }

      purrr::pmap(sdGp, function(resName, windowSdName, topVarNum, topVarSamples, topVarNumInput, topVarSamplesInput, inputChanged) {

        if (!is.na(windowSdName)) {
          th <- dataTable %>%
            dplyr::pull({{windowSdName}}) %>%
            dplyr::nth(topVarNum)

          dataTable <- dataTable %>%
            dplyr::filter(!!dplyr::sym(windowSdName) >= th)

          if (length(topVarSamples) == length(samples)) {

            if(verbose) {
              message(glue::glue("------------------------------
                             Filtering windows based on standard deviation across all {length(topVarSamples)} {groupString}s ({windowSdName}).
                             Standard deviation threshold = {format(th, digits = 3)} resulting in {nrow(dataTable)} windows."))
            }
          } else {

            if(verbose) {
              message(glue::glue("------------------------------
                             Filtering windows based on standard deviation across {length(topVarSamples)} {groupString}s ({windowSdName}).
                             Standard deviation threshold = {format(th, digits = 3)} resulting in {nrow(dataTable)} windows."))
            }
          }

        } else {

          if(verbose){
            message(glue::glue("------------------------------
                             No filtering of windows based on window standard deviation."))
          }
          th <- NA
        }

        dataTable <- dataTable %>%
          dplyr::mutate(window = getWindowNames(.)) %>%
          tibble::column_to_rownames("window") %>%
          dplyr::select(tidyr::all_of(samples))

        if(verbose) {
          message(glue::glue("Performing {method} with {ncol(dataTable)} {groupString}s and {nrow(dataTable)} windows
                           -> {resName}."))
        }

        dataTable <- dataTable %>%
          t()

        if (!is.numeric(dataTable)) {
          stop("Input contains non-numeric values.")
        }

        if (any(is.na(dataTable))) {
          stop("Input contains missing values.")
        }

        if (any(!is.finite(dataTable))) {
          stop("Input contains infinite values.")
        }

        if (method == "UMAP") {

          uwotObj <- dataTable %>%
            uwot::umap(n_components = 3,
                       ...) %>%
            as.data.frame()

          colnames(uwotObj) <- paste0("UMAP",seq_along(1:ncol(uwotObj)))

          return(list(resObj =  list(x = uwotObj, windows = colnames(dataTable)), th = th))

        } else if (method == "PCA") {

          prcompObj <- dataTable %>%
            stats::prcomp(center = center, scale. = scale, rank. = nPC)

          return(list(resObj = prcompObj, th = th))

        } else {
          stop(glue::glue("Method {method} not known! Options are PCA or UMAP."))
        }

      }) %>%
        purrr::set_names(sdGp$resName)
    }) %>%
    purrr::list_flatten()

  if(verbose) { message("------------------------------") }

  th <- purrr::map(res, "th") %>%
    unlist()

  res <- purrr::map(res, "resObj")

  if (method == "UMAP") {
    windows <- purrr::map(res, "windows")
  } else if (method == "PCA") {
    windows <- purrr::map(res, ~ rownames(.x$rotation))
  }

  paramList <- list(method = method,
                    regionsToOverlap = regionsToOverlap,
                    normMethod = normMethod,
                    minEnrichment = minEnrichment,
                    useGroupMeans = useGroupMeans,
                    minDensity = minDensity,
                    topVar = topVar,
                    windowSdThreshold = th)

  if (method == "UMAP") {
    paramList <- c(paramList, ...)

    elements <- purrr::map2(res, windows, function(x,y) {
      methods::new("mesaUMAP",
          points = x$x,
          windows = y
      )
    }
    ) %>%
      rlang::set_names(nm = names(res))

  } else if (method == "PCA") {
    paramList <- c(paramList, list(center = center,
                                   scale = scale,
                                   nPC = nPC))

    elements <- purrr::map2(res, windows, function(x,y) {
      methods::new("mesaPCA",
          prcomp = x,
          windows = y
      )
    }
    ) %>%
      rlang::set_names(nm = names(res))

  }

  windowFilteringList <- list(initial = initialNumWindows,
                              notInRegionsToOverlap = numWindowsRemovedRegionOverlap,
                              belowMinDensity = numWindowsRemovedMinDensity,
                              containMissingVals = numWindowsRemovedMissingVals)

  out <- methods::new("mesaDimRed",
              samples = samples,
              sampleTable = qseaSet %>% qsea::getSampleTable(),
              params = c(paramList,windowFilteringList),
              dataTable = if (returnDataTable) { dataTable } else { data.frame() },
              res = elements
  )

  return(out)
}


getShapeScale <- function(plotData, shape, shapePalette = NULL, colourScaleType = NULL, NAshape = NULL) {
  
  if (shape == "NULLshape") {
    return(ggplot2::scale_shape_manual(values = shapePalette, guide = "none"))
  }
  
  shapeVals <- dplyr::pull(plotData, {{shape}})
  nShape    <- shapeVals %>% setdiff(NA) %>% unique() %>% length()
  hasNA     <- any(is.na(shapeVals))
  
  if (is.null(shapePalette)) {
    if (!is.null(colourScaleType) && colourScaleType == "diverging") {
      shapePalette <- if (hasNA) c(21, 24, 22, 23) else c(21, 24, 22, 23, 25)
      if (nShape > length(shapePalette)) {
        stop(glue::glue("`shape` variable '{shape}' has {nShape} unique values; the maximum allowed by default when using a divergent colour scale is {length(shapePalette)} unique values."))
      }

    } else if (nShape <= 4 && hasNA) {
      shapePalette <- c(21, 24, 22, 23)
    } else if (nShape <= 5 && !hasNA) {
      shapePalette <- c(21, 24, 22, 23, 25)
    } else if (hasNA) {
      shapePalette <- c(16, 4, 0, 17, 8, 9, 15, 13, 2, 18, 14, 3, 1, 5, 6, 10, 11, 12)
    } else {
      shapePalette <- c(16, 4, 0, 17, 8, 9, 15, 13, 2, 18, 14, 3, 1, 5, 6, 10, 11, 12, 7)
    }
    if (nShape > length(shapePalette)) {
      stop(glue::glue("`shape` variable '{shape}' has {nShape} unique values; {if (hasNA) 'with' else 'without'} NAs present, the maximum allowed by default is {length(shapePalette)}."))
    }
  } else if (is.numeric(shapePalette)) {
    if (nShape > length(shapePalette)) {
      stop(glue::glue("`shape` variable '{shape}' has {nShape} unique values; `shapePalette` only has {length(shapePalette)} values."))
    }
    if (any(0:20 %in% shapePalette) & any(21:25 %in% shapePalette)) {
      stop("'shapePalette' must contain either 0-20 (unfilled) OR 21-25 (filled), not both.")
    }
    
  } else if (is.character(shapePalette)) {
    shapesInput <- shapePalette
    if (shapePalette == "filled+border") {
      shapePalette <- if (hasNA) c(21, 24, 22, 23) else c(21, 24, 22, 23, 25)
    } else {
      shapePalette <- switch(
        shapePalette,
        "line-first"   = c(1, 8, 2, 0, 9, 3, 13, 6, 14, 4, 5, 10, 11, 12, 16, 17, 15, 18),
        "filled-first" = c(16, 17, 15, 18, 1, 8, 2, 0, 9, 3, 13, 6, 14, 4, 5, 10, 11, 12),
        "mixture"      = c(16, 4, 0, 17, 8, 9, 15, 13, 2, 18, 14, 3, 1, 5, 6, 10, 11, 12),
        stop("`shapePalette` must be 'line-first', 'filled-first', 'mixture', 'filled+border', numeric, or NULL.")
      )
      if (colourScaleType == "diverging") {
        warning(glue::glue("Using '{shapesInput}' palette with a divergent colour scale may reduce visibility near zero."))
      }
    }
    if (nShape > length(shapePalette)) {
      stop(glue::glue("`shape` variable '{shape}' has {nShape} unique values; only {length(shapePalette)} shapes available with `shapePalette` = '{shapesInput}'."))
    }
    
  } else {
    stop("`shapePalette` must be 'line-first', 'filled-first', 'mixture', 'filled+border', numeric, or NULL.")
  }
  
  if (hasNA) {
    usingFilled <- any(21:25 %in% shapePalette)
    defaultNA   <- if (usingFilled) 25 else 7
    defaultInUse <- defaultNA %in% shapePalette[1:nShape]
    
    if (is.null(NAshape)) {
      if (defaultInUse) {
        stop(glue::glue(
          "The default NAshape = {defaultNA} is already in use within the shapePalette. ",
          "Specify a different NAshape or remove the default shape from the shapePalette."
        ))
      }
      NAshape <- defaultNA
      
    } else if ((!usingFilled && NAshape %in% 21:25) || (usingFilled && NAshape %in% 0:20)) {
      warning(glue::glue(
        "'shapePalette' is using {if (usingFilled) 'filled' else 'unfilled'} shapes, ",
        "but specified NAshape = {NAshape} is from the {if (!usingFilled) 'filled' else 'unfilled'} set. ",
        "Resetting NAshape to default = {defaultNA} for {if (usingFilled) 'filled' else 'unfilled'} shapes."
      ))
      if (defaultInUse) {
        stop(glue::glue("The default NAshape = {defaultNA} is also already in use within the shapePalette."))
      }
      NAshape <- defaultNA
      
    } else if (NAshape %in% shapePalette[1:nShape]) {
      paletteName <- if (exists("shapesInput")) shapesInput else NULL
      
      if (!is.null(paletteName) && paletteName %in% c("line-first", "filled-first")) {
        if (defaultInUse) {
          stop(glue::glue(
            "Both the specified NAshape = {NAshape} and the default NAshape = {defaultNA} ",
            "are already in use within the shapePalette."
          ))
        }
        warning(glue::glue(
          "Specified NAshape = {NAshape} is already in use in the shapePalette. ",
          "Because palette = '{paletteName}' has a fixed order, NAshape has been reset to default = {defaultNA}."
        ))
        NAshape <- defaultNA
      } else {
        if (defaultInUse) {
          stop(glue::glue(
            "Both the specified NAshape = {NAshape} and the default NAshape = {defaultNA} ",
            "are already in use within the shapePalette."
          ))
        }
        warning(glue::glue(
          "Specified NAshape = {NAshape} is already in use in the shapePalette. ",
          "Replacing that entry in the shapePalette with the default NAshape = {defaultNA} at the end, ",
          "and keeping NAshape = {NAshape} for missing values."
        ))
        clash_idx <- match(NAshape, shapePalette)
        shapePalette <- c(shapePalette[-clash_idx], defaultNA)
      }
    }
    
  } else {
    NAshape <- NA
  }
  
  
  ggplot2::scale_shape_manual(values = shapePalette, na.value = NAshape)
}


getGeomPoint <- function(cV, shape, my_scale_shape, pointSize, alpha) {

  if(any(my_scale_shape$palette(1) %in% 21:25)){
    filledShapes <- TRUE
  } else {
    filledShapes <- FALSE 
  }
  
  if (filledShapes) {
    my_geom_point <- ggplot2::geom_point(ggplot2::aes(fill = !!rlang::sym(cV), shape = !!rlang::sym(shape)),
                                         colour = "black", size = pointSize, alpha = alpha)
  } else {
    my_geom_point <- ggplot2::geom_point(ggplot2::aes(colour = !!rlang::sym(cV), shape = !!rlang::sym(shape)),
                                         size = pointSize, alpha = alpha)
  }

  return(my_geom_point)

}

getColourScale <- function(plotData, cV, cols, colourScaleType, my_scale_shape, NAcolour, symDivColourScale) {

  if(any(my_scale_shape$palette(1) %in% 21:25)){
    filledShapes <- TRUE
  } else {
    filledShapes <- FALSE 
  }
  
  if (length(cV) == 1 && cV == "NULLcol") {
    my_scale_colour <- if (filledShapes) {
      ggplot2::scale_fill_identity()
    } else {
      ggplot2::scale_colour_identity()
    }
  } else {
    if (is.null(cols)) {
      if (colourScaleType == "qualitative") {
        my_scale_colour <- if (filledShapes) {
          hues::scale_fill_iwanthue(na.value = NAcolour)
        } else {
          hues::scale_colour_iwanthue(na.value = NAcolour)
        }
      } else if (colourScaleType == "sequential_non_neg") {
        my_scale_colour <- if (filledShapes) {
          ggplot2::scale_fill_viridis_c(na.value = NAcolour)
        } else {
          ggplot2::scale_colour_viridis_c(na.value = NAcolour)
        }
      } else if (colourScaleType == "sequential_non_pos") {
        my_scale_colour <- if (filledShapes) {
          ggplot2::scale_fill_viridis_c(direction = -1, na.value = NAcolour)
        } else {
          ggplot2::scale_colour_viridis_c(direction = -1, na.value = NAcolour)
        }
      } else if (colourScaleType == "diverging") {
        cols <- RColorBrewer::brewer.pal(9, "RdBu") %>% rev()
        cols[5] <- "grey90"
      }
      
    } else {
      if (colourScaleType == "qualitative") {
        nCol <- plotData %>% pull(cV) %>% setdiff(NA) %>% unique() %>% length()
        if (nCol > length(cols)) {
          stop(glue::glue("`colour` variable '{cV}' has {nCol} unique values; the `colourPalette` argument only has {length(cols)} unique values."))
        }
        my_scale_colour <- if (filledShapes) {
          ggplot2::scale_fill_manual(values = cols, na.value = NAcolour)
        } else {
          ggplot2::scale_colour_manual(values = cols, na.value = NAcolour)
        }
        
      } else if (colourScaleType == "sequential_non_neg" | colourScaleType == "sequential_non_pos") {
        my_scale_colour <- if (filledShapes) {
          ggplot2::scale_fill_gradientn(colours = cols, na.value = NAcolour)
        } else {
          ggplot2::scale_colour_gradientn(colours = cols, na.value = NAcolour)
        }
      }
    }
    
    if (colourScaleType == "diverging") {
      
      cVdat <- plotData[[cV]]
      
      if (symDivColourScale) {
        maxCV <- max(cVdat, na.rm = TRUE)
        minCV <- min(cVdat, na.rm = TRUE)
        absMinCV <- abs(minCV)
        if (abs(minCV) < maxCV) {
          # minimum (negative) value is smaller in magnitude than the largest (positive) value; colour scale needs to be extended beyond the minimum value
          vals <- scales::rescale(c(-maxCV, 0, maxCV),
                                  to = c(-(maxCV - absMinCV) / (absMinCV + maxCV), 1))
        } else {
          # minimum (negative) value is larger in magnitude than the largest (positive) value; colour scale needs to be extended beyond the maximum value
          vals <- scales::rescale(c(minCV, 0, -minCV),
                                  to = c(0, 1 +  (absMinCV - maxCV) / (absMinCV + maxCV)))
        }
        
        my_scale_colour <- if (filledShapes) {
          ggplot2::scale_fill_gradientn(colours = cols, values = vals, na.value = NAcolour)
        } else {
          ggplot2::scale_colour_gradientn(colours = cols, values = vals, na.value = NAcolour)
        }
        
      } else {
        my_scale_colour <- if (filledShapes) {
          ggplot2::scale_fill_gradientn(colours = cols,
                                        values = scales::rescale(c(min(cVdat, na.rm = TRUE), 0, max(cVdat, na.rm = TRUE))),
                                        na.value = NAcolour)
        } else {
          ggplot2::scale_colour_gradientn(colours = cols,
                                          values = scales::rescale(c(min(cVdat, na.rm = TRUE), 0, max(cVdat, na.rm = TRUE))),
                                          na.value = NAcolour)
        }
      }
    }
  }

  return(my_scale_colour)

}

getLegendParams <- function(cV, shape, my_scale_shape, colourScaleType) {
  filledShapes <- any(my_scale_shape$palette(1) %in% 21:25)
  
  if (filledShapes &&
      !(cV %in% c("NULLcol", shape)) &&
      colourScaleType == "qualitative") {
    my_legend_params <- ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21, col = "black")))
  } else {
    my_legend_params <- NULL
  }
  return(my_legend_params)
}

#' Plot principal component analysis (PCA) results
#'
#' Generate publication-ready plots from the output of [getPCA()]. The function
#' supports flexible visual encodings (colour, shape, size) of sample-level
#' metadata and can return multiple plots across different PC pairs and
#' annotation variables.
#'
#' @param object `mesaDimRed`.  
#'   A dimensionality-reduction container returned by [getPCA()].
#' @param components `integer(2)` **or** `list` of `integer(2)`.  
#'   PC indices to plot (e.g. `c(1, 2)` for PC1 vs PC2), or a list of such pairs
#'   to produce multiple plots.  
#'   **Default:** `list(c(1, 2), c(2, 3))`.
#'
#' @param colour `NULL` or `character()` .  
#'   Name(s) of sample-table variable(s) mapped to point colour; one plot per
#'   variable.  
#'   **Default:** `NULL` (no colour mapping).
#'
#' @param colourPalette `NULL` or `character()`.  
#'   Palette used for colouring points. If `NULL`, a suitable default is chosen.  
#'   **Default:** `NULL`.
#'
#' @param NAcolour `character(1)`.  
#'   Colour for missing values in `colour`.  
#'   **Default:** `"grey50"`.
#'
#' @param symDivColourScale `logical(1)`.  
#'   If `TRUE`, diverging colour scales are centred symmetrically around zero.
#'   Ignored for non-diverging scales.  
#'   **Default:** `FALSE`.
#'
#' @param shape `NULL` or `character(1)`.  
#'   Sample-table variable mapped to point shape (single variable only).  
#'   **Default:** `NULL`.
#'   
#' @param shapePalette `NULL`, `character(1)`, or `numeric()`.  
#'   Shapes used for categories of the `shape` variable. Options:  
#'   * A numeric vector of base‐R shape codes.  
#'     - 0–20: line/filled shapes.  
#'     - 21–25: filled shapes with a border (border colour = black).  
#'   * A keyword, with internally defined sets:  
#'     - `"line-first"`: 15 line shapes, then 4 filled shapes (**max 19 categories**).  
#'     - `"filled-first"`: 4 filled shapes, then 15 line shapes (**max 19 categories**).  
#'     - `"mixture"`: mixture of line and filled shapes (**max 19 categories**).  
#'     - `"filled+border"`: filled+border shapes (**max 5 categories, or 4 if NAs in `shape`**).  
#'   * `NULL`: choose automatically (`"mixture"` for non-diverging/none; `"filled+border"`
#'     for diverging colour scales).  
#'   **Default:** `NULL`.
#'   
#' @param NAshape `NULL` or `numeric(1)`.  
#'   Shape for missing values in `shape`. If `NULL`, uses `7`, or `25` when
#'   `"filled+border"` is active.  
#'   **Default:** `NULL`.
#'
#' @param showSampleNames `logical(1)`.  
#'   Overlay sample names on points.  
#'   **Default:** `FALSE`.
#'
#' @param pointSize `numeric(1)`.  
#'   Point size.  
#'   **Default:** `2`.
#'
#' @param alpha `numeric(1)`.  
#'   Point transparency in `[0, 1]`.  
#'   **Default:** `1`.
#'
#' @param plotlyAnnotations `character()` .  
#'   Column names from the sample table to add as tooltips when converting to
#'   interactive `plotly`.  
#'   **Default:** `""` (empty string).
#'
#' @return A `list` of `ggplot` objects—one per combination of `components`
#'   pairs and `colour` variables—each depicting the chosen PCs with requested
#'   aesthetics. Specifically:
#'   * one plot per `components` pair × per `colour` variable;
#'   * shapes/labels applied per `shape`/`showSampleNames`;
#'   * consistent axis labelling (PC variance if available).
#'
#' @details
#' **Shape defaults**  
#' * `"mixture"` is used for non-diverging or discrete colour scales (or no colour).  
#' * `"filled+border"` is used for diverging colour scales.  
#' These defaults aim to keep categories visually separable for larger cohorts.
#'
#' @seealso [getPCA()], [plotUMAP()], [plotDimRed()], [qsea::getSampleTable()]
#' @family dimred-plotting
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Basic PCA plot (PC1 vs PC2) with defaults
#' exampleTumourNormal %>%
#'   getPCA() %>%
#'   plotPCA()
#'
#' # Colour by sample group; show names
#' exampleTumourNormal %>%
#'   getPCA() %>%
#'   plotPCA(colour = "group",
#'           showSampleNames = TRUE)
#'
#' # Plot PC1 vs PC3; shape by tissue
#' exampleTumourNormal %>%
#'   getPCA() %>%
#'   plotPCA(components = list(c(1, 3)),
#'           colour = "group",
#'           shape  = "tissue") 
#'           
#' @export
plotPCA.mesaDimRed <- function(object,
                    components = list(c(1, 2), c(2, 3)),
                    colour = NULL,
                    colourPalette = NULL,
                    NAcolour = "grey50",
                    symDivColourScale = FALSE,
                    shape = NULL,
                    shapePalette = NULL,
                    NAshape = NULL,
                    showSampleNames = FALSE,
                    pointSize = 2,
                    alpha = 1,
                    plotlyAnnotations = "") {

  out <- plotDimRed(object = object,
             components = components,
             colour = {{colour}},
             colourPalette = colourPalette,
             NAcolour = NAcolour,
             symDivColourScale = symDivColourScale,
             shape = {{shape}},
             shapePalette = shapePalette,
             NAshape = NAshape,
             showSampleNames = showSampleNames,
             pointSize = pointSize,
             alpha = alpha,
             plotlyAnnotations = plotlyAnnotations)

  return(out)
}

#' Plot UMAP results
#'
#' Convenience wrapper around [plotDimRed()] for UMAP embeddings returned by
#' [getUMAP()]. Supports flexible visual encodings (colour, shape, size) of
#' sample-level metadata, and can return multiple plots across UMAP component
#' pairs and annotation variables.
#'
#' @param object `mesaDimRed`.  
#'   A dimensionality-reduction container returned by [getUMAP()].
#'
#' @param components `integer(2)` **or** `list` of `integer(2)`.  
#'   UMAP indices to plot (e.g. `c(1, 2)` for UMAP1 vs UMAP2), or a list of such
#'   pairs to produce multiple plots.  
#'   **Default:** `list(c(1, 2))`.
#'
#' @param colour `NULL` or `character()`.  
#'   Name(s) of sample-table variable(s) mapped to point colour; one plot per
#'   variable.  
#'   **Default:** `NULL` (no colour mapping).
#'
#' @param colourPalette `NULL` or `character()`.  
#'   Palette used for colouring points. If `NULL`, a suitable default is chosen.  
#'   **Default:** `NULL`.
#'
#' @param NAcolour `character(1)`.  
#'   Colour for missing values in `colour`.  
#'   **Default:** `"grey50"`.
#'
#' @param symDivColourScale `logical(1)`.  
#'   If `TRUE`, diverging colour scales are centred symmetrically around zero.
#'   Ignored for non-diverging scales.  
#'   **Default:** `FALSE`.
#'
#' @param shape `NULL` or `character(1)`.  
#'   Sample-table variable mapped to point shape (single variable only).  
#'   **Default:** `NULL`.
#'
#' @param shapePalette `NULL`, `character(1)`, or `numeric()`.  
#'   Shapes used for categories of the `shape` variable. Options:  
#'   * A numeric vector of base-R shape codes.  
#'     - 0–20: line/filled shapes.  
#'     - 21–25: filled shapes with a border (border colour = black).  
#'   * A keyword, with internally defined sets:  
#'     - `"line-first"`: 15 line shapes, then 4 filled shapes (**max 19 categories**).  
#'     - `"filled-first"`: 4 filled shapes, then 15 line shapes (**max 19 categories**).  
#'     - `"mixture"`: mixture of line and filled shapes (**max 19 categories**).  
#'     - `"filled+border"`: filled+border shapes (**max 5 categories, or 4 if NAs in `shape`**).  
#'   * `NULL`: choose automatically (`"mixture"` for non-diverging/none; `"filled+border"`
#'     for diverging colour scales).  
#'   **Default:** `NULL`.
#'
#' @param NAshape `NULL` or `numeric(1)`.  
#'   Shape for missing values in `shape`. If `NULL`, uses `7`, or `25` when
#'   `"filled+border"` is active.  
#'   **Default:** `NULL`.
#'
#' @param showSampleNames `logical(1)`.  
#'   Overlay sample names on points.  
#'   **Default:** `FALSE`.
#'
#' @param pointSize `numeric(1)`.  
#'   Point size.  
#'   **Default:** `2`.
#'
#' @param alpha `numeric(1)`.  
#'   Point transparency in `[0, 1]`.  
#'   **Default:** `1`.
#'
#' @param plotlyAnnotations `character()`.  
#'   Column names from the sample table to add as tooltips when converting to
#'   interactive `plotly`.  
#'   **Default:** `""` (empty string).
#'
#' @return A `list` of `ggplot` objects—one per combination of `components`
#'   pairs and `colour` variables—each depicting the chosen UMAP dimensions with
#'   requested aesthetics. Specifically:
#'   * one plot per `components` pair × per `colour` variable;  
#'   * shapes/labels applied per `shape`/`showSampleNames`;  
#'   * consistent axis labelling (UMAP1/UMAP2 etc.).
#'
#' @family dimred-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Default UMAP (beta-normalised, top 1000 most variable windows)
#' exampleTumourNormal %>%
#'   getUMAP(n_neighbors = 5, min_dist = 1) %>%
#'   plotUMAP()
#'
#' # Colour by group; show names
#' exampleTumourNormal %>%
#'   getUMAP(n_neighbors = 5, min_dist = 1) %>%
#'   plotUMAP(colour = "group",
#'            showSampleNames = TRUE)
#'
#' # Custom components and shape mapping
#' exampleTumourNormal %>%
#'   getUMAP(n_neighbors = 5, min_dist = 1) %>%
#'   plotUMAP(components = list(c(1, 2)),
#'            colour = "group",
#'            shape  = "tissue")
#'
#' @export
plotUMAP <- function(object,
                     components = list(c(1, 2)),
                     colour = NULL,
                     colourPalette = NULL,
                     NAcolour = "grey50",
                     symDivColourScale = FALSE,
                     shape = NULL,
                     shapePalette = NULL,
                     NAshape = NULL,
                     showSampleNames = FALSE,
                     pointSize = 2,
                     alpha = 1,
                     plotlyAnnotations = "") {

  out <- plotDimRed(object = object,
             components = components,
             colour = colour,
             colourPalette = colourPalette,
             NAcolour = NAcolour,
             symDivColourScale = symDivColourScale,
             shape = shape,
             shapePalette = shapePalette,
             NAshape = NAshape,
             showSampleNames = showSampleNames,
             pointSize = pointSize,
             alpha = alpha,
             plotlyAnnotations = plotlyAnnotations)

  return(out)
}


#' Plot dimensionality reduction results
#'
#' Visualise PCA or UMAP coordinates stored in a [mesaDimRed] object.  
#' Wrappers such as [plotPCA()] and [plotUMAP()] provide convenient shortcuts
#' for common use cases.
#'
#' @param object `mesaDimRed`.  
#'   A dimensionality reduction container returned by [getPCA()] or [getUMAP()].
#'
#' @param components `integer(2)` **or** `list` of `integer(2)`.  
#'   Component indices to plot (e.g. `c(1, 2)` for PC1 vs PC2). A list of pairs
#'   generates multiple plots.  
#'   **Default:** `list(c(1, 2), c(2, 3))` for PCA, `list(c(1, 2))` for UMAP.
#'
#' @param colour `NULL` or `character()`.  
#'   Name(s) of sample-table variable(s) mapped to point colour. One plot is
#'   produced per variable.  
#'   **Default:** `NULL` (no colouring).
#'
#' @param colourPalette `NULL` or `character()`.  
#'   Palette for point colours. If `NULL`, defaults are selected automatically.  
#'   **Default:** `NULL`.
#'
#' @param NAcolour `character(1)`.  
#'   Colour for missing values in `colour`.  
#'   **Default:** `"grey50"`.
#'
#' @param symDivColourScale `logical(1)`.  
#'   If `TRUE`, diverging colour scales are centred symmetrically around zero.  
#'   Ignored for non-diverging scales.  
#'   **Default:** `FALSE`.
#'
#' @param shape `NULL` or `character(1)`.  
#'   Sample-table variable mapped to point shape. Only one variable supported.  
#'   **Default:** `NULL`.
#'
#' @param shapePalette `NULL`, `character(1)`, or `numeric()`.  
#'   Shapes used for categories of the `shape` variable. Options:  
#'   * A numeric vector of base-R shape codes:  
#'     - 0–20: line/filled shapes.  
#'     - 21–25: filled shapes with borders (border colour = black).  
#'   * A keyword:  
#'     - `"line-first"`: 15 line shapes, then 4 filled shapes (**max 19 categories**).  
#'     - `"filled-first"`: 4 filled shapes, then 15 line shapes (**max 19 categories**).  
#'     - `"mixture"`: mixture of line and filled shapes (**max 19 categories**).  
#'     - `"filled+border"`: filled+border shapes (**max 5 categories; 4 if NAs in `shape`**).  
#'   * `NULL`: automatic choice—`"mixture"` unless a diverging colour scale is
#'     used, in which case `"filled+border"`.  
#'   **Default:** `NULL`.
#'
#' @param NAshape `NULL` or `numeric(1)`.  
#'   Shape used for missing values in `shape`. If `NULL`, defaults to `7`, or `25`
#'   when `"filled+border"` is active.  
#'   **Default:** `NULL`.
#'
#' @param showSampleNames `logical(1)`.  
#'   If `TRUE`, overlay sample names on points.  
#'   **Default:** `FALSE`.
#'
#' @param pointSize `numeric(1)`.  
#'   Size of plotted points.  
#'   **Default:** `2`.
#'
#' @param alpha `numeric(1)`.  
#'   Transparency of points in `[0, 1]`.  
#'   **Default:** `1`.
#'
#' @param plotlyAnnotations `character()`.  
#'   Column names from the sample table used as tooltips when converting plots
#'   to interactive `plotly`.  
#'   **Default:** `""` (empty string).
#'
#' @return A `list` of `ggplot2` objects, one for each combination of:  
#'   * component pairs in `components`, and  
#'   * variables specified in `colour`.  
#'
#' Each plot depicts the requested dimensions with aesthetics mapped as specified.
#'
#' @family dimred-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # PCA: colour by group
#' exampleTumourNormal %>%
#'   getPCA(nPC = 3) %>%
#'   plotDimRed(colour = "group")
#'
#' # UMAP: colour by group, shape by tissue
#' exampleTumourNormal %>%
#'   getUMAP(n_neighbors = 5, min_dist = 1) %>%
#'   plotDimRed(colour = "group",
#'              shape = "tissue") 
#'
#' # Show sample names
#' exampleTumourNormal %>%
#'   getPCA(nPC = 3) %>%
#'   plotDimRed(colour = "group",
#'              showSampleNames = TRUE) 
#'
#' @export
plotDimRed <- function(object,
                    components = list(c(1, 2), c(2, 3)),
                    colour = NULL,
                    colourPalette = NULL,
                    NAcolour = "grey50",
                    symDivColourScale = FALSE,
                    shape = NULL,
                    shapePalette = NULL,
                    NAshape = NULL,
                    showSampleNames = FALSE,
                    pointSize = 2,
                    alpha = 1,
                    plotlyAnnotations = ""
){

  if (!inherits(object,"mesaDimRed")) {
    stop("First argument should be the output from the getPCA()/getUMAP() functions in mesa.")
  }
  
  sampleTable <- object@sampleTable
  
  if (!is.null(colour)){
    colDiff <- setdiff(colour, colnames(sampleTable))
    if(length(colDiff) > 0){
      stop(glue::glue("Can't colour by {colDiff}, as it is not present in the sampleTable!
                      
                      ")) #empty line is required here
    }
  }

  if (!is.null(shape)){
    colDiff <- setdiff(shape, colnames(sampleTable))
    if(length(colDiff) > 0){
      stop(glue::glue("Can't set shapes by {colDiff}, as it is not present in the sampleTable!
                      
                      ")) #empty line is required here
    }
  }
    
  if (!is.list(components)) {
    components <- list(components)
  }

  if (object@params$method == "PCA") {
    columnPrefix <- "PC"
  } else if (object@params$method == "UMAP") {
    columnPrefix <- "UMAP"
  } else {
    stop("Method {object@params$method} not known")
  }

  if (length(plotlyAnnotations) > 1) {
    plotlyAnnotations <- plotlyAnnotations %>% purrr::set_names(., nm = .)
  } else if (plotlyAnnotations != "") {
    plotlyAnnotations <- plotlyAnnotations %>% purrr::set_names(., nm = .)
  }

  components <- components %>% purrr::set_names(purrr::map(components, ~ glue::glue("{columnPrefix}{.x}") %>% glue::glue_collapse("vs")))

  ggp <- purrr::imap(object@res, function(single, resName) {

    numWindows <- length(single@windows)

    if (columnPrefix == "PC") {
      propVar <- single@prcomp$sdev ^ 2 / sum(single@prcomp$sdev ^ 2)
      propVar <- round(propVar * 100, 2)
      plotData <- single@prcomp$x
    } else {
      propVar = NULL
      plotData <- single@points
    }

    plotData <- plotData %>%
      tibble::as_tibble(rownames = "sample_name") %>%
      dplyr::left_join(sampleTable, by = "sample_name")

    if (!is.null(colourPalette) & is.null(colour)) {
      stop("`colourPalette` argument is non-NULL, but `colour` argument is NULL.")
    }

    if (!is.null(shape) & length(shape) > 1) {
      stop("Argument `shape` can only be of length one.")
    }

    if (is.null(shape)) {
      plotData <- plotData %>%
        dplyr::mutate(NULLshape = "21")

      if (is.null(shapePalette)){
        shapePalette = 21
      }
      
      shape <- "NULLshape"
    }

    if (is.null(colour)) {
      plotData <- plotData %>%
        dplyr::mutate(NULLcol = "black")

      colour <- "NULLcol"
    }

    makePlot <- function(components, plotData, numWindows, my_geom_point, my_scale_colour, my_scale_shape, my_legend_params) {

      env <- new.env(parent = globalenv())
      env$plotData <- plotData
      env$columnPrefix <- columnPrefix
      env$components <- components
      env$plotlyAnnotations <- plotlyAnnotations

      topVarInfo <- object@params$topVar %>% filter(resName == !!resName)

      if (is.na(topVarInfo$topVarNum)) {
        titleString <- glue::glue("all {numWindows} windows")
        subtitleString <- glue::glue("Using {object@params$normMethod} values.")
      } else {
        titleString <- glue::glue("top {numWindows} most variable windows")
        if (length(topVarInfo$topVarSamples[[1]]) == length(object@samples)) {
          titleSubstring <- "all "
        } else {
          titleSubstring <- ""
        }
        subtitleString <- glue::glue("Using {object@params$normMethod} values and {titleSubstring}{length(topVarInfo$topVarSamples[[1]])} samples to calculate std dev.")
      }

      ggp <- with(env, {
        ggplot2::ggplot(plotData,
                        ggplot2::aes(!!rlang::sym(glue::glue("{columnPrefix}{components[1]}")),
                                     !!rlang::sym(glue::glue("{columnPrefix}{components[2]}")),
                                     label = sample_name,
                                     !!!rlang::syms(plotlyAnnotations)
                        ))
      }) +
        my_geom_point +
        my_scale_colour +
        my_scale_shape +
        my_legend_params +
        ggplot2::ggtitle(glue::glue("{object@params$method} for {length(object@samples)} samples using {titleString}."),
                         subtitle = subtitleString) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 12.5))

      if (object@params$method == "PCA") {
        ggp <- ggp +
          ggplot2::xlab(glue::glue("PC{components[1]} ({propVar[components[1]]}%)")) +
          ggplot2::ylab(glue::glue("PC{components[2]} ({propVar[components[2]]}%)"))
      }

      if (showSampleNames) {
        if (!requireNamespace("ggrepel", quietly = TRUE)) {
          message("Package \"ggrepel\" is recommended to repel labels. Using default method.")
          ggp <- ggp + ggplot2::geom_text()
        } else {
          ggp <- ggp + ggrepel::geom_text_repel()
        }
      }

      return(ggp)
    }

    if (length(colour) == 1 && colour == "NULLcol") {

      my_scale_shape <- getShapeScale(plotData, shape, shapePalette, NAshape = NAshape)

      my_geom_point <- getGeomPoint(colour, shape, my_scale_shape, pointSize = pointSize, alpha = alpha)

      my_scale_colour <- getColourScale(cV = colour, my_scale_shape = my_scale_shape)
      
      my_legend_params <- getLegendParams(colour, shape, my_scale_shape)

      ggp <- purrr::map(components, makePlot, plotData, numWindows, my_geom_point, my_scale_colour, my_scale_shape, my_legend_params)

      return(ggp)

    } else {

      ggp <- purrr::map2(purrr::set_names(colour), list(colourPalette), function(cV, cols) {

        cVdat <- plotData[[cV]]

        colourScaleType <- dplyr::case_when(is.factor(cVdat) || is.character(cVdat) ~ "qualitative", # qualitative variable
                                            is.numeric(cVdat) && min(cVdat, na.rm = TRUE) >= 0 ~ "sequential_non_neg", # non-negative sequential variable
                                            is.numeric(cVdat) && max(cVdat, na.rm = TRUE) <= 0 ~ "sequential_non_pos", # non-positive sequential variable
                                            is.numeric(cVdat) && (max(cVdat, na.rm = TRUE) > 0 && min(cVdat, na.rm = TRUE) < 0) ~ "diverging") # diverging variable

        if (is.na(colourScaleType)) {
          stop(glue::glue("The variable `{cV}` can not be mapped to a colour scale."))
        }

        my_scale_shape <- getShapeScale(plotData, shape, shapePalette, colourScaleType, NAshape = NAshape)

        my_geom_point <- getGeomPoint(cV, shape, my_scale_shape, pointSize = pointSize, alpha = alpha)

        my_scale_colour <- getColourScale(plotData, cV, cols, colourScaleType, my_scale_shape, NAcolour = NAcolour, symDivColourScale = symDivColourScale)
        
        my_legend_params <- getLegendParams(cV, shape, my_scale_shape, colourScaleType)

        ggp <- purrr::map(components, makePlot, plotData, numWindows, my_geom_point, my_scale_colour, my_scale_shape, my_legend_params)

        return(ggp)

      })
    }
  })

  return(ggp)
}
