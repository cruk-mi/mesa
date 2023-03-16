#' This function takes the output of [getPCA()] and produces PCA plots.
#' @param object The output from [getPCA()].
#' @param qseaSet The qseaSet object used to generate `object`.
#' @param components Vector of the two components to plot, or a list of vectors to make multiple plots. Default is to produce plots for PC1 vs PC2 and PC2 vs PC3.
#' @param colour Character vector of variable names from the qseaSet sample table for setting the colour of the points (samples). Separate plots are made for each variable.
#' @param colourPalette Character vector giving the colour palette to use for the points (samples). Defaults are used if not supplied.
#' @param NAcolour Colour to use for NA values in the `colour` variable. Default is "grey50".
#' @param symDivColourScale Logical indicating if a diverging colour scale should be symmetric around zero (ignored if a diverging colour scale is not used).
#' @param shape Character giving variable name from the qseaSet sample table for setting the shape of the points (samples). Can only accept a single variable name.
#' @param shapePalette Shapes to use for the points (samples) for each category of the `shape` variable. Can be one of the following options:
#' * A numeric vector specifying the set of shapes. Can be either integers between 0 and 20 (line or filled shapes) or integers between 21 and 25 (filled shapes with a border; border colour is set to black).
#' * A character specifying which types of shapes to use (with the exact set of shapes set internally by the function). Either:
#'     * "line-first" (15 line shapes, then 4 filled shapes; max. 19 categories).
#'     * "filled-first" (4 filled shapes, then 15 line shapes; max. 19 categories).
#'     * "mixture" (mixture of line and filled shapes; max. 19 categories).
#'     * "filled+border" (max. 5 categories, or 4 if there are NAs in the `shape` variable).
#' * NULL; defaults are used which is the "mixture" set of shapes for non-diverging colour scales (or no colour scale) and "filled+border" for diverging colour scales.
#' @param NAshape Shape to use for NA values in the `shape` variable. Default is shape 7, or shape 25 if filled shapes with a border are being used.
#' @param showSampleNames Logical indicating whether to show the sample names.
#' @param pointSize Numeric value to set the size of the points
#' @return A ggplot object or list of ggplot objects
#' @export
#'
plotPCA <- function(object,
                    qseaSet,
                    components = list(c(1, 2)),
                    colour = NULL,
                    colourPalette = NULL,
                    NAcolour = "grey50",
                    symDivColourScale = FALSE,
                    shape = NULL,
                    shapePalette = NULL,
                    NAshape = NULL,
                    showSampleNames = FALSE,
                    pointSize = 2
){

  if (!("pcas" %in% names(object))) {
    stop("First argument should be a pca object from getPCA")
  }

  if (!is.qseaSet(qseaSet)) {
    stop("Second argument should be a qseaSet")
  }

  if (!is.list(components)) {
    components <- list(components)
  }

  components <- components %>% purrr::set_names(purrr::map(components, ~ glue::glue("PC{.x}") %>% glue::glue_collapse("vs")))

  ggp <- purrr::imap(object$pcas, function(pca, pcaName) {

    propVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    propVar <- round(propVar * 100, 2)

    plotData <- pca$x %>%
      tibble::as_tibble(rownames = "sample_name") %>%
      dplyr::left_join(qsea::getSampleTable(qseaSet))

    if (!is.null(colourPalette) & is.null(colour)) {
      stop("`colourPalette` argument is non-NULL, but `colour` argument is NULL.")
    }

    if (!is.null(shapePalette) & is.null(shape)) {
      stop("`shapePalette` argument is non-NULL, but `shape` argument is NULL.")
    }

    if (!is.null(shape) & length(shape) > 1) {
      stop("Argument `shape` can only be of length one.")
    }

    if (is.null(shape)) {
      plotData <- plotData %>%
        dplyr::mutate(NULLshape = 16)

      shape <- "NULLshape"
    }

    if (is.null(colour)) {
      plotData <- plotData %>%
        dplyr::mutate(NULLcol = "black")

      colour <- "NULLcol"
    }

    getShapeScale <- function(plotData, shape, shapePalette, colourScaleType = NULL) {

      if (is.null(NAshape)) {
        NAshape <- 7
      }

      if (shape == "NULLshape") {
        my_scale_shape <- ggplot2::scale_shape_identity(na.value = NAshape)

      } else {

        nShape <- plotData %>% pull(shape) %>% setdiff(NA) %>%  unique() %>% length()

        if (is.null(shapePalette)) {
          if (!is.null(colourScaleType) && colourScaleType == "diverging") {
            if (any(is.na(plotData %>% pull(shape)))) {
              shapePalette <- c(21, 24, 22, 23)
              NAshape <- 25
            } else {
              shapePalette <- c(21, 24, 22, 23, 25)
            }
            if (nShape > length(shapePalette)) {
              stop(glue::glue("`shape` variable '{shape}' has {nShape} unique values; the maximum allowed by default when using a divergent colour scale is {length(shapePalette)} unique values."))
            }

          } else {
            shapePalette <- c(16, 8, 0, 17, 3, 9, 15, 13, 2, 18, 14, 4, 1, 5, 6, 10, 11, 12)
            if (nShape > length(shapePalette)) {
              stop(glue::glue("`shape` variable '{shape}' has {nShape} unique values; the maximum allowed by default is {length(shapePalette)} unique values."))
            }
          }

        } else if (is.numeric(shapePalette)) {
          if (nShape > length(shapePalette)) {
            stop(glue::glue("`shape` variable '{shape}' has {nShape} unique values; the `shapePalette` argument only has {length(shapePalette)} unique values."))
          }
          if (any(0:20 %in% shapePalette) & any(21:25 %in% shapePalette)) {
            stop("'shapePalette' argument can either contain integers between 0 and 20 (line or filled shapes), or between 21 and 25 (filled shapes with borders), but not both.")
          }
          if (any(21:25 %in% shapePalette)) {
            NAshape <- 25
          }

        } else if (is.character(shapePalette)) {
          shapesInput <- shapePalette
          if (shapePalette == "filled+border") {
            if (any(is.na(plotData %>% pull(shape)))) {
              shapePalette <- c(21, 24, 22, 23)
              NAshape <- 25
            } else {
              shapePalette <- c(21, 24, 22, 23, 25)
            }
          } else {
            if (shapePalette == "line-first") {
              shapePalette <- c(1, 8, 2, 0, 9, 3, 13, 6, 14, 4, 5, 10, 11, 12, 16, 17, 15, 18)
            } else if (shapePalette == "filled-first") {
              shapePalette <- c(16, 17, 15, 18, 1, 8, 2, 0, 9, 3, 13, 6, 14, 4, 5, 10, 11, 12)
            } else if (shapePalette == "mixture") {
              shapePalette <- c(16, 8, 0, 17, 3, 9, 15, 13, 2, 18, 14, 4, 1, 5, 6, 10, 11, 12)
            } else {
              stop("`shapePalette` argument can take the following character values: 'line-first', 'filled-first', 'mixture' or 'filled+border'; or can be numeric or NULL.")
            }

            if (colourScaleType == "diverging") {
              warning(glue::glue("Using the '{shapesInput}' colour scale with a divergent colour scale may lead to points around zero on the colour scale being almost invisible."))
            }
          }
          if (nShape > length(shapePalette)) {
            stop(glue::glue("`shape` variable '{shape}' has {nShape} unique values; there are only {length(shapePalette)} shapes available with argument `shapePalette` = '{shapesInput}'."))
          }
        } else {
          stop("`shapePalette` argument is not in a valid format. It can take the following character values: 'line-first', 'filled-first', 'mixture' or 'filled+border'; or can be numeric or NULL.")
        }

        if (NAshape %in% shapePalette[1:nShape] & any(is.na(plotData %>% pull(shape)))) {
          stop(glue::glue("NA shape value (={NAshape}) is already being used for a '{shape}' category. Values in use: {paste0(shapePalette[1:nShape], collapse = ', ')}."))
        }

        my_scale_shape <- ggplot2::scale_shape_manual(values = shapePalette, na.value = NAshape)
      }

      return(my_scale_shape)

    }

    getGeomPoint <- function(cV, shape, my_scale_shape) {

      filledShapes <- ifelse(my_scale_shape$scale_name == "manual" & any(my_scale_shape$palette(1) %in% 21:25),
                             TRUE, FALSE)

      if (filledShapes) {
        my_geom_point <- ggplot2::geom_point(ggplot2::aes(fill = !!rlang::sym(cV), shape = !!rlang::sym(shape)),
                                             colour = "black", size = pointSize)
      } else {
        my_geom_point <- ggplot2::geom_point(ggplot2::aes(colour = !!rlang::sym(cV), shape = !!rlang::sym(shape)),
                                             size = pointSize)
      }

      return(my_geom_point)

    }

    getColourScale <- function(plotData, cV, cols, colourScaleType, my_scale_shape) {

      filledShapes <- ifelse(my_scale_shape$scale_name == "manual" & any(my_scale_shape$palette(1) %in% 21:25),
                             TRUE, FALSE)

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

      return(my_scale_colour)

    }

    makePlot <- function(PCs, plotData, my_geom_point, my_scale_colour, my_scale_shape) {

      env <- new.env(parent = globalenv())
      env$plotData <- plotData
      env$PCs <- PCs

      topVarInfo <- object$params$topVar %>% filter(pcaName == !!pcaName)

      if (is.na(topVarInfo$topVarNum)) {
        titleString <- glue::glue("all {length(object$windows[[pcaName]])} windows")
        subtitleString <- glue::glue("Using {object$params$normMethod} values.")
      } else {
        titleString <- glue::glue("top {length(object$windows[[pcaName]])} most variable windows")
        if (length(topVarInfo$topVarSamples[[1]]) == length(object$samples)) {
          titleSubstring <- "all "
        } else {
          titleSubstring <- ""
        }
        subtitleString <- glue::glue("Using {object$params$normMethod} values and {titleSubstring}{length(topVarInfo$topVarSamples[[1]])} samples to calculate std dev.")
      }

      ggp <- with(env, {
        ggplot2::ggplot(plotData,
                        ggplot2::aes(!!rlang::sym(glue::glue("PC{PCs[1]}")), !!rlang::sym(glue::glue("PC{PCs[2]}")),
                                     label = sample_name))
      }) +
        my_geom_point +
        my_scale_colour +
        my_scale_shape +
        ggplot2::xlab(glue::glue("PC{PCs[1]} ({propVar[PCs[1]]}%)")) +
        ggplot2::ylab(glue::glue("PC{PCs[2]} ({propVar[PCs[2]]}%)")) +
        ggplot2::ggtitle(glue::glue("PCA for {length(object$samples)} samples using {titleString}."),
                         subtitle = subtitleString) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 12.5))

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

      my_scale_shape <- getShapeScale(plotData, shape, shapePalette)

      my_geom_point <- getGeomPoint(colour, shape, my_scale_shape)

      my_scale_colour <- ggplot2::scale_colour_identity()

      ggp <- purrr::map(components, makePlot, plotData, my_geom_point, my_scale_colour, my_scale_shape)

      return(ggp)

    } else {

      ggp <- purrr::map2(purrr::set_names(colour), list(colourPalette), function(cV, cols) {

        cVdat <- plotData[[cV]]

        colourScaleType <- dplyr::case_when(is.factor(cVdat) | is.character(cVdat) ~ "qualitative", # qualitative variable
                                            is.numeric(cVdat) & min(cVdat, na.rm = TRUE) >= 0 ~ "sequential_non_neg", # non-negative sequential variable
                                            is.numeric(cVdat) & max(cVdat, na.rm = TRUE) <= 0 ~ "sequential_non_pos", # non-positive sequential variable
                                            is.numeric(cVdat) & (max(cVdat, na.rm = TRUE) > 0 & min(cVdat, na.rm = TRUE) < 0) ~ "diverging") # diverging variable

        if (is.na(colourScaleType)) {
          stop(glue::glue("The variable `{cV}` can not be mapped to a colour scale."))
        }

        my_scale_shape <- getShapeScale(plotData, shape, shapePalette, colourScaleType)

        my_geom_point <- getGeomPoint(cV, shape, my_scale_shape)

        my_scale_colour <- getColourScale(plotData, cV, cols, colourScaleType, my_scale_shape)

        ggp <- purrr::map(components, makePlot, plotData, my_geom_point, my_scale_colour, my_scale_shape)

        return(ggp)

      })
    }
  })

  return(ggp)
}


#' This function is a modified version of the [qsea::getPCA()] function
#' @param qseaSet A qseaSet object.
#' @param dataTable A data frame of normalised values for a set of windows (rows) and samples (columns), e.g. from [getDataTable()]. It must have seqnames, start and end columns. Can also be a [GenomicRanges::GRanges()] object with normalised values in the metadata columns.
#' @param regionsToOverlap Optional. Only windows in x overlapping `regionsToOverlap` will be considered. A [GenomicRanges::GRanges()] object or a data frame which can be coerced into a [GenomicRanges::GRanges()] object.
#' @param normMethod What normalisation method to use. Typically a character giving name of predefined normalisation method (e.g. "beta" or "nrpm"). See [qsea::normMethod()].
#' @param minEnrichment Minimum number of reads for beta values to not give NA. Passed to [getDataTable()].
#' @param useGroupMeans Whether to average samples over the group column (i.e. combine replicates)
#' @param minDensity A minimum CpG density level to filter out windows with values lower than.
#' @param topVarNum Number of most variable windows to keep. Defaults to 1000. If value is NA, NULL, Inf or is at least as big as the available number of windows, then all available windows are used for the PCA. Can also be a vector, in which case PCA is performed for each component.
#' @param topVarSamples Samples to use to determine variability. Either NULL or NA (use all samples; default), a character vector of sample names or a regular expression to match sample names on. Can also be a list, in which case PCA is performed for each list component. If `length(topVarNum) > 1`, list should be same length as `topVarNum` and each component will be matched with corresponding component of `topVarNum`.
#' @param center A logical value indicating if windows should be centred to have mean zero. Default is TRUE.
#' @param scale A logical value indicating if windows should be scaled to have unit variance. Default is FALSE.
#' @param nPC Number of principal components to be calculated. Default is 5.
#' @param returnDataTable A logical value indicating if the table of normalised values (generated by [getDataTable()] or by processing the input `dataTable`) should be returned as part of the output. See below for details.
#' @return A list with the following components:
#' \item{pcas}{A named list of objects of class prcomp. List names correspond to entries in the `params$topVar` component of the output.}
#' \item{samples}{Names of the samples analysed.}
#' \item{windows}{A list of names of the windows analysed. List names correspond to entries in the `params$topVar` component of the output.}
#' \item{dataTable}{If `returnDataTable = TRUE`, the table of normalised values. The table will contain windows after filtering based on `regionsToOverlap` and `minDensity`, but prior to filtering based on `topVarNum`. Windows with missing values will also have been removed. The table will also include column(s) of window standard deviations, where applicable, with column name(s) corresponding to entries in the `params$topVar` component of the output.}
#' \item{params}{A list of the parameters used.}
#' \item{windowFiltering}{A list of the initial number of windows and the number filtered out at each filtering step (not including filtering based on `topVarNum`).}
#' @export
#'
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
                   returnDataTable = FALSE) {

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
      warning("Newly-generated dataTable from qseaSet on a small subset of samples/windows does not match input dataTable.
              \n")
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

  message(glue::glue("------------------------------
                     Initial number of windows = {initialNumWindows}."))

  if (!is.null(regionsToOverlap)) {
    regionsToOverlap <- regionsToOverlap %>%
      tibble::as_tibble() %>%
      dplyr::select(tidyselect::any_of(c("seqnames", "start", "end", "CpG_density"))) %>%  # keep only minimum columns necessary
      asValidGranges()

    if (is.null(dataTable)) {
      qseaSet <- qseaSet %>%
        filterByOverlaps(regionsToOverlap = regionsToOverlap)

      numWindowsRemovedRegionOverlap <- initialNumWindows - length(getWindows(qseaSet))
      message(glue::glue("Filtered out {numWindowsRemovedRegionOverlap} windows using regionsToOverlap: {length(getRegions(qseaSet))} windows remaining."))

    } else {
      dataTable <- dataTable %>%
        plyranges::filter_by_overlaps(y = regionsToOverlap)

      numWindowsRemovedRegionOverlap <- initialNumWindows - length(dataTable)
      message(glue::glue("Filtered out {numWindowsRemovedRegionOverlap} windows using regionsToOverlap: {length(dataTable)} windows remaining."))
    }
  } else {
    numWindowsRemovedRegionOverlap <- NULL
  }

  if (is.null(dataTable)) {
    message("-----------")
    dataTable <- qseaSet %>%
      filterWindows(CpG_density >= minDensity) %>%
      getDataTable(normMethod = normMethod, useGroupMeans = useGroupMeans)
    message("-----------")

    if (minDensity > 0) {
      numWindowsRemovedMinDensity <- length(getWindows(qseaSet)) - nrow(dataTable)
      message(glue::glue("Filtered out {numWindowsRemovedMinDensity} windows with CpG_density < {minDensity}: {nrow(dataTable)} windows remaining."))
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
                         dplyr::select(seqnames, start, end, CpG_density)) %>%
      dplyr::filter(CpG_density >= minDensity)

    if (minDensity > 0) {
      numWindowsRemovedMinDensity <- currentNumWindows - nrow(dataTable)
      message(glue::glue("Filtered out {numWindowsRemovedMinDensity} windows with CpG_density < {minDensity}: {nrow(dataTable)} windows remaining."))
    } else {
      numWindowsRemovedMinDensity <- NULL
    }

  }

  currentNumWindows <- nrow(dataTable)

  dataTable <- dataTable %>%
    tidyr::drop_na(tidyr::all_of(samples))

  numWindowsRemovedMissingVals <- currentNumWindows - nrow(dataTable)

  if (numWindowsRemovedMissingVals > 0) {
    message(glue::glue("Filtered out {numWindowsRemovedMissingVals} windows with at least one missing value: {nrow(dataTable)} windows remaining.\n
                       ------------------------------
                       \n"))
  } else {
    message(glue::glue("No windows have missing values.\n
            ------------------------------
            \n"))
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
      message("These values are not used; PCA is already being done with all remaining windows.\n")
    } else {
      message("These values are not used; PCA will be done with all remaining windows instead.\n")
    }
  }


  topVar <- topVar %>%
    dplyr::distinct(topVarNum, topVarSamples, .keep_all = TRUE) %>%
    dplyr::group_by(topVarSamples) %>%
    dplyr::mutate(windowSdName = ifelse(!is.na(topVarNum), glue::glue("windowSd{dplyr::cur_group_id()}"), NA), .before = 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(windowSdName, topVarNum) %>%
    dplyr::mutate(pcaName = glue::glue("pca{dplyr::row_number()}"), .before = 1) %>%
    dplyr::bind_rows(topVar) %>%
    dplyr::distinct(topVarNum, topVarSamples, topVarNumInput, topVarSamplesInput, .keep_all = TRUE)

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
          message(glue::glue("Calculating standard deviation for each window across all {length(samples)} {groupString}s:
                           {paste0(tVS, collapse = ', ')}.
                           -> column name {nm}.
                           \n"))
        } else {
          message(glue::glue("Calculating standard deviation for each window across {length(tVS)} of {length(samples)} {groupString}s:
                           {paste0(tVS, collapse = ', ')}.
                           -> column name {nm}.
                           \n"))
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
      dplyr::left_join(dataTable, .)
  }

  pca <- topVar %>%
    dplyr::filter(!is.na(pcaName)) %>%
    dplyr::group_by(windowSdName) %>%
    dplyr::group_map(.keep = TRUE, .f = function(sdGp, sdName) {

      if (!is.na(sdName$windowSdName)) {
        dataTable <- dataTable %>%
          dplyr::arrange(dplyr::desc(!!dplyr::sym(sdName$windowSdName)))
      }

      purrr::pmap(sdGp, function(pcaName, windowSdName, topVarNum, topVarSamples, ...) {

        if (!is.na(windowSdName)) {
          th <- dataTable %>%
            dplyr::pull({{windowSdName}}) %>%
            dplyr::nth(topVarNum)

          dataTable <- dataTable %>%
            dplyr::filter(!!dplyr::sym(windowSdName) >= th)

          if (length(topVarSamples) == length(samples)) {

            message(glue::glue("------------------------------\n
                             Filtering windows based on standard deviation across all {length(topVarSamples)} {groupString}s ({windowSdName}).
                             Standard deviation threshold = {format(th, digits = 3)} resulting in {nrow(dataTable)} windows.
                             \n"))
          } else {
            message(glue::glue("------------------------------\n
                             Filtering windows based on standard deviation across {length(topVarSamples)} {groupString}s ({windowSdName}).
                             Standard deviation threshold = {format(th, digits = 3)} resulting in {nrow(dataTable)} windows.
                             \n"))
          }


        } else {

          message(glue::glue("------------------------------\n
                             No filtering of windows based on window standard deviation.
                             \n"))
          th <- NA
        }

        dataTable <- dataTable %>%
          dplyr::mutate(window = getWindowNames(.)) %>%
          tibble::column_to_rownames("window") %>%
          dplyr::select(tidyr::all_of(samples))

        message(glue::glue("Performing PCA with {ncol(dataTable)} {groupString}s and {nrow(dataTable)} windows
                           -> {pcaName}.
                             \n"))

        dataTable <- dataTable %>%
          t()

        if (!is.numeric(dataTable)) {
          stop("Input to PCA contains non-numeric values.")
        }

        if (any(is.na(dataTable))) {
          stop("Input to PCA contains missing values.")
        }

        if (any(!is.finite(dataTable))) {
          stop("Input to PCA contains infinite values.s")
        }

        prcompObj <- dataTable %>%
          stats::prcomp(center = center, scale. = scale, rank. = nPC)

        return(list(prcompObj = prcompObj, th = th))

      }) %>%
        purrr::set_names(sdGp$pcaName)
    }) %>%
    purrr::list_flatten()

  message("------------------------------")

  th <- purrr::map(pca, "th") %>%
    unlist()

  pcas <- purrr::map(pca, "prcompObj")

  windows <- purrr::map(pcas, ~ rownames(.x$rotation))

  return(list(pcas = pcas,
              samples = samples,
              windows = windows,
              dataTable = if (returnDataTable) {
                dataTable
              } else {
                NULL
              },
              params = list(regionsToOverlap = regionsToOverlap,
                            normMethod = normMethod,
                            minEnrichment = minEnrichment,
                            useGroupMeans = useGroupMeans,
                            minDensity = minDensity,
                            topVar = topVar,
                            windowSdThreshold = th,
                            center = center,
                            scale = scale,
                            nPC = nPC),
              windowFiltering = list(initial = initialNumWindows,
                                     notInRegionsToOverlap = numWindowsRemovedRegionOverlap,
                                     belowMinDensity = numWindowsRemovedMinDensity,
                                     containMissingVals = numWindowsRemovedMissingVals)))
}