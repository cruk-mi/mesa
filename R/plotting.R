#' Heatmap of signal across selected genomic regions
#'
#' Plot sample signal over a set of regions as a heatmap using **ComplexHeatmap**.
#' Regions are selected by overlap with `regionsToOverlap`, values are taken from
#' the `qseaSet`, and optional annotations from the sample/region metadata can be added.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing windows and counts/betas.
#'
#' @param regionsToOverlap `GRanges`, `data.frame` (coercible to `GRanges`), or `NULL`.  
#'   If provided, only windows overlapping these regions are plotted; otherwise all windows are considered.  
#'   **Default:** `NULL`.
#'
#' @param normMethod `character(1)`.  
#'   Which measure to plot, `"beta"` or `"nrpm"`.  
#'   **Default:** `"beta"`.
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, average samples by the `group` column in the sample table (combine replicates).  
#'   **Default:** `FALSE`.
#'
#' @param sampleAnnotation `character()` or `NULL`.  
#'   Names of columns from the sample table to display as column annotations (e.g., `c("tumour","tissue")`).  
#'   **Default:** `NULL`.
#'
#' @param windowAnnotation `character()` or `NULL`.  
#'   Names of columns from region metadata (from the `qseaSet` regions or `regionsToOverlap`) to display as row annotations (e.g., `c("CpG_density")`).  
#'   **Default:** `NULL`.
#'
#' @param clusterNum `integer(1)` or `NULL`.  
#'   If set, cut the **column** dendrogram into this many clusters and display cluster labels.  
#'   **Default:** `NULL`.
#'
#' @param maxScale `numeric(1)`.  
#'   Upper limit of the colour scale (used for `"nrpm"`; not applied to `"beta"`).  
#'   **Default:** `5`.
#'
#' @param clip `numeric(1)`.  
#'   Cap values above this threshold before plotting (applies to `"nrpm"`, ignored for `"beta"`).  
#'   **Default:** `1e9`.
#'
#' @param minDensity `numeric(1)`.  
#'   Minimum `CpG_density` required to keep a window.  
#'   **Default:** `0`.
#'
#' @param clusterRows `logical(1)`.  
#'   Whether to cluster rows (windows).  
#'   **Default:** `FALSE`.
#'
#' @param clusterCols `logical(1)`.  
#'   Whether to cluster columns (samples).  
#'   **Default:** `TRUE`.
#'
#' @param clusterMethod `character(1)`.  
#'   Clustering method for dendrograms (e.g., `"ward.D2"`).  
#'   **Default:** `"ward.D2"`.
#'
#' @param annotationColors `list` or `NA`.  
#'   Optional named colour maps for annotations, e.g. `list(tumour = c(Tumour = "firebrick4", Normal = "blue"))`.  
#'   **Default:** `NA`.
#'
#' @param minEnrichment `numeric(1)`.  
#'   For `"beta"`, values with enrichment `< minEnrichment` are set to `NA`.  
#'   **Default:** `3`.
#'
#' @param showSampleNames `logical(1)` or `NULL`.  
#'   If `NULL`, names are shown when there are fewer than 50 samples; set `TRUE`/`FALSE` to force.  
#'   **Default:** `NULL`.
#'
#' @param annotationPosition `character(1)`.  
#'   Where to place column annotations (e.g., `"right"` or `"bottom"`).  
#'   **Default:** `"right"`.
#'
#' @param title `character(1)` or `NULL`.  
#'   Optional plot title.  
#'   **Default:** `NULL`.
#'
#' @param ... Additional arguments forwarded to **ComplexHeatmap** constructors.
#'
#' @details
#' Windows are filtered by overlap (if `regionsToOverlap` is given) and by `minDensity`.
#' For `"beta"`, low-enrichment values are set to `NA` using `minEnrichment`. For `"nrpm"`,
#' `clip` and `maxScale` help stabilise colour scaling for outliers. If `useGroupMeans = TRUE`,
#' samples are aggregated by `group` before plotting.
#'
#' @return
#' Draws a heatmap on the active device and (invisibly) returns the underlying
#' **ComplexHeatmap** object (e.g., a `Heatmap`/`HeatmapList`) for further composition with
#' [ComplexHeatmap::draw()].
#'
#' @seealso
#' [qsea::makeTable()], [ComplexHeatmap::Heatmap()], [ComplexHeatmap::draw()],
#' [calculateDMRs()], [getSampleTable()], [getWindows()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' 
#' # Compute regions (DMRs) to plot
#' DMRs <- exampleTumourNormal %>% calculateDMRs(variable = "tumour", contrasts = "first")
#' 
#' # Basic heatmap of beta values over DMRs
#' exampleTumourNormal %>% plotRegionsHeatmap(DMRs)
#' 
#' # Cluster rows, add sample and window annotations
#' exampleTumourNormal %>% 
#'   plotRegionsHeatmap(DMRs, 
#'                      clusterRows = TRUE, 
#'                      sampleAnnotation = c(tumour, tissue))
#' 
#' # Group means, 2 column clusters, and custom annotation colours
#' exampleTumourNormal %>% 
#'   plotRegionsHeatmap(regionsToOverlap = DMRs, 
#'                    clusterRows = TRUE, 
#'                    clusterNum = 2,
#'                    sampleAnnotation = tumour,
#'                    windowAnnotation = CpG_density,
#'                    annotationColors = list(tumour = c("Tumour" = "firebrick4", "Normal" = "blue"))
#'                     )
#' 
#' @export
#' @examples
#' # calculate DMRs to plot
#' DMRs <- exampleTumourNormal %>% calculateDMRs(variable = "tumour", contrasts = "first")
#' # plot these windows
#' exampleTumourNormal %>% plotRegionsHeatmap(DMRs)
#' # cluster the rows and add annotation
#' exampleTumourNormal %>% plotRegionsHeatmap(DMRs, clusterRows = TRUE, sampleAnnotation = c(tumour, tissue))
#' # more complex example
#' exampleTumourNormal %>% 
#'   plotRegionsHeatmap(regionsToOverlap = DMRs, 
#'                    clusterRows = TRUE, 
#'                    clusterNum = 2,
#'                    sampleAnnotation = tumour,
#'                    windowAnnotation = CpG_density,
#'                    annotationColors = list(tumour = c("Tumour" = "firebrick4", "Normal" = "blue"))
#'                     )
#' 
plotRegionsHeatmap <- function(qseaSet, regionsToOverlap = NULL,
                                normMethod = "beta",
                                sampleAnnotation = NULL,
                                windowAnnotation = NULL,
                                annotationColors = NA,
                                useGroupMeans = FALSE,
                                clusterRows = FALSE,
                                clusterCols = TRUE,
                                minEnrichment = 3,
                                maxScale = 5,
                                clusterNum = NULL,
                                clip = 1000000000,
                                minDensity = 0,
                                annotationPosition = "right",
                                title = NULL,
                                showSampleNames = NULL,
                                clusterMethod = "ward.D2", ...) {

  if (is.null(regionsToOverlap)) {
    regionsToOverlap <- qseaSet %>%
      qsea::getRegions()
  }

  regionsToOverlap <- asValidGranges(regionsToOverlap)

  if (length(regionsToOverlap) > 20000) {stop("More than 20000 regions requested.")}

  if (length(regionsToOverlap) == 0) {stop("No genomic regions given!")}

  if (normMethod == "beta") {maxScale = min(clip,1)}

  col_fun = circlize::colorRamp2(seq(0, maxScale, length.out = 9), RColorBrewer::brewer.pal(name = "YlOrRd", n = 9))

  clipFn <- function(x, a, b) {a + (x - a > 0) * (x - a) - (x - b > 0) * (x - b)}

  if (!is.list(annotationColors)) {

    namesVec <- names(annotationColors)
    names(namesVec) <- namesVec

    purrr::map(namesVec, function(x){
      usedValues <- qseaSet %>%
        qsea::getSampleTable() %>%
        dplyr::pull(x) %>%
        unique()
      return(x = annotationColors[[x]][usedValues])
    }
    )

  }

  # define a function that removes rows that have 1 row.
  remove_almost_empty_rows <- function(dat)  {
    mask_keep <- rowSums(is.na(dat)) != (ncol(dat) - 1)
    janitor:::remove_message(dat = dat, mask_keep = mask_keep, which = "rows", reason = "almost empty")
    return(dat[mask_keep, , drop = FALSE])
  }

  dataTab <- qseaSet %>%
    filterByOverlaps(regionsToOverlap = regionsToOverlap) %>%
    filterWindows(CpG_density >= minDensity) %>%
    getDataTable(normMethod = normMethod,
                 useGroupMeans = useGroupMeans,
                 minEnrichment = minEnrichment
    ) %>%
    dplyr::mutate(window = paste0(seqnames, ":",start, "-",end)) %>%
    dplyr::mutate_all( ~ dplyr::case_when(!is.nan(.x) ~ .x)) # do something with NaN values?

  if (useGroupMeans) {
    colsToFind <- qseaSet %>% getSampleGroups2() %>% names()
  } else {
    colsToFind <- qseaSet %>% qsea::getSampleNames()
  }

  numData <- dataTab %>%
    tibble::column_to_rownames("window") %>%
    dplyr::select(tidyselect::all_of(colsToFind)) %>%
    clipFn(a = 0, b = clip) %>%
    janitor::remove_empty(which = "cols", quiet = FALSE) %>%
    janitor::remove_empty(which = "rows", quiet = FALSE)

  if (clusterRows) {
    numData <- numData %>%
      remove_almost_empty_rows()
  }

  dataTab <- dataTab %>%
    filter(window %in% rownames(numData))

  windowAnnotationDf <- getWindowAnnotation(dataTab, regions = regionsToOverlap,
                                            windowAnnotation = {{windowAnnotation}},
                                            clusterRows = clusterRows)

  # Ensure regions order is maintained if not clustering.
  numData <- numData[rownames(windowAnnotationDf),]

  if (!useGroupMeans) {

    annots <- qseaSet %>%
      filter(sample_name %in% !!(colnames(numData))) %>%
      makeHeatmapAnnotations(useGroupMeans = useGroupMeans,
                             specifiedAnnotationColors = annotationColors,
                             sampleAnnotation = {{sampleAnnotation}},
                             windowAnnotationDf = windowAnnotationDf)

  } else {
    annots <- qseaSet %>%
      filter(group %in% !!(colnames(numData))) %>%
      makeHeatmapAnnotations(useGroupMeans = useGroupMeans,
                             specifiedAnnotationColors = annotationColors,
                             sampleAnnotation = {{sampleAnnotation}},
                             windowAnnotationDf = windowAnnotationDf )
  }

  colAnnot <- annots$sample
  rowAnnot <- annots$window

  if (ncol(dataTab) == 1) {
    clusterCols <- FALSE
  }

  if (clusterRows) {
    dataTab <- remove_almost_empty_rows(dataTab)
  }

  if (clusterCols && !is.null(clusterNum) && clusterNum > 1) {
    colSplit <- clusterNum
  } else {
    colSplit <- NULL
  }

  annotName <- dplyr::case_when(normMethod == "beta" ~ "Beta value",
                                normMethod == "nrpm" ~ "NRPM",
                                TRUE ~ stringr::str_to_title(normMethod)
  )

  if (is.null(showSampleNames)) {
    if (ncol(numData) > 50) {
      showSampleNames <- FALSE
    }  else {
      showSampleNames <- TRUE
    }
  }

  if(clusterRows){
    if(numData %>% stats::dist() %>% is.na() %>% sum() > 0) {
      message("Can not cluster rows due to too many missing values, setting clusterRows to be FALSE")
      clusterRows <- FALSE
    }
  }

  numData %>%
    as.matrix() %>%
    ComplexHeatmap::Heatmap(name = annotName,
                            cluster_rows = clusterRows,
                            cluster_columns = clusterCols,
                            show_row_names = FALSE,
                            col = col_fun,
                            clustering_method_rows = clusterMethod,
                            clustering_method_columns = clusterMethod,
                            heatmap_legend_param = list(legend_direction = "horizontal",
                                                        at = seq(0, maxScale, length.out = 6) %>% round(1)),
                            column_split = colSplit,
                            show_column_names = showSampleNames,
                            column_title = NULL,
                            top_annotation = colAnnot,
                            left_annotation = rowAnnot,
                            ...) %>%
    ComplexHeatmap::draw(heatmap_legend_side = "bottom",
                         annotation_legend_side = annotationPosition,
                         column_title = title,
                         column_title_gp = grid::gpar(fontsize = 16))

}


#' Prepare window-level annotations for plotRegionsHeatmap()
#'
#' Note that windowAnnotation must be enclosed within double curly brackets when used.
#'
#' Build a data frame of row annotations aligned to the windows shown in the
#' heatmap created by [plotRegionsHeatmap()]. Annotations can be sourced from
#' the table used to plot (e.g. `getDataTable()` output) and/or from metadata
#' columns on a `GRanges` of regions.
#'
#' @param dataTab `data.frame`/`tibble`.  
#'   Table produced by [getDataTable()] that includes `seqnames`, `start`, `end`,
#'   and a `window` column of the form `"chr:start-end"`.
#'
#' @param regions `GRanges`.  
#'   Genomic regions used for overlap-derived annotations (e.g., columns already
#'   present in `mcols(regions)` will be joinable by window coordinates).
#'
#' @param windowAnnotation tidyselect specification or `character()` or `NULL`.  
#'   Columns to keep as annotations. Supports tidyselect helpers when supplied
#'   unquoted (e.g., `CpG_density`, `starts_with("qc_")`) or a character vector
#'   of column names.  
#'   **Default:** `NULL` (keep none).
#'
#' @param clusterRows `logical(1)`.  
#'   If `FALSE`, preserve genomic order of windows. If `TRUE`, follow the row
#'   order produced by clustering in the heatmap.  
#'   **Default:** `FALSE`.
#'
#' @return `data.frame` with row names set to the `window` identifiers and the
#'   requested annotation columns.
#'
#' @details
#' Rows are aligned by the `window` label (`"chr:start-end"`). When both `dataTab`
#' and `regions` provide overlapping information, columns from either source may be
#' selected via `windowAnnotation`. Use tidyselect to pick or pattern-match columns.
#'
#' @seealso
#' [plotRegionsHeatmap()], [makeHeatmapAnnotations()], [getDataTable()]
#'
#' @keywords internal
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' exampleTumourNormal %>%
#'   filterByOverlaps(qsea::getRegions(exampleTumourNormal)[1:50]) %>%
#'   getDataTable(normMethod = "nrpm") %>%
#'   dplyr::mutate(window = paste0(seqnames, ":", start, "-", end)) %>%
#'   mesa:::getWindowAnnotation(
#'     regions = qsea::getRegions(exampleTumourNormal)[1:50],
#'     windowAnnotation = CpG_density
#'   ) %>%
#'   head()
#'   
getWindowAnnotation <- function(dataTab, regions, windowAnnotation = NULL, clusterRows = FALSE) {
  rowAnnotDf <- dataTab %>%
    plyranges::as_granges() %>%
    plyranges::join_overlap_left(regions %>%
                                   tibble::as_tibble() %>%
                                   tibble::rowid_to_column(".rowID") %>%
                                   plyranges::as_granges(),
                                 suffix = c("",".regionsToOverlap")) %>%
    tibble::as_tibble() %>%
    dplyr::select(window, .rowID, {{windowAnnotation}})

  if(!clusterRows) {
    rowAnnotDf <- rowAnnotDf  %>%
      dplyr::arrange(.rowID)
  }

  rowAnnotDfdistinct <- rowAnnotDf %>%
    dplyr::select(-.rowID) %>%
    dplyr::distinct()
  
  
  nonUniqueWindowsDf <- rowAnnotDfdistinct %>% 
    dplyr::group_by(window) %>% 
    dplyr::filter(dplyr::n() > 1)
  
  if(nrow(nonUniqueWindowsDf) > 0) {
    stop(glue::glue("Non-unique annotations found in window annotations: 
         {paste(capture.output(print(head(nonUniqueWindowsDf))), collapse = '\n')}
                    "))
  }
  rowAnnotDfdistinct %>%
    tibble::column_to_rownames("window")

}

#' Create ComplexHeatmap annotation objects (samples and windows)
#'
#' Helper to generate `ComplexHeatmap::HeatmapAnnotation()` objects for samples
#' and windows with automatic colour handling for categorical and numeric variables.
#'
#' @param qseaSet `qseaSet`.  
#'   The input object from which sample metadata are taken.
#'
#' @param sampleAnnotation tidyselect specification or `character()` or `NULL`.  
#'   Columns from the sample table to annotate (unquoted helpers like `group`, `tumour`,
#'   or a character vector of names).  
#'   **Default:** `NULL`.
#'
#' @param windowAnnotationDf `data.frame` or `NULL`.  
#'   Row annotations aligned to windows (e.g., output of [getWindowAnnotation()]).  
#'   **Default:** `NULL`.
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, annotate groups instead of individual samples (using the `group`
#'   column of the sample table).  
#'   **Default:** `FALSE`.
#'
#' @param specifiedAnnotationColors `list` or `NA`.  
#'   Optional mapping of levels → colours to override auto-generated palettes,
#'   e.g. `list(tumour = c(Tumour = "firebrick4", Normal = "blue"))`.  
#'   **Default:** `NA`.
#'
#' @param windowOrientation `character(1)`.  
#'   Orientation for window (row) annotations: `"row"` or `"column"`.  
#'   **Default:** `"row"`.
#'
#' @param sampleOrientation `character(1)`.  
#'   Orientation for sample (column) annotations: `"column"` or `"row"`.  
#'   **Default:** `"column"`.
#'
#' @return `list(sample = <HeatmapAnnotation|NULL>, window = <HeatmapAnnotation|NULL>)`
#'   ready to pass to `ComplexHeatmap::Heatmap(...)`.
#'
#' @details
#' When `specifiedAnnotationColors` is not supplied, discrete variables receive a
#' qualitative palette and numeric variables are mapped to a continuous gradient.
#' The function respects `useGroupMeans` to annotate either per-sample or per-group.
#'
#' @seealso
#' [plotRegionsHeatmap()], [getWindowAnnotation()], [ComplexHeatmap::HeatmapAnnotation()],
#' [ComplexHeatmap::Heatmap()]
#'
#' @keywords internal
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Build window annotations and immediately feed them into makeHeatmapAnnotations()
#' exampleTumourNormal %>%
#'   filterByOverlaps(qsea::getRegions(exampleTumourNormal)[1:50]) %>%
#'   getDataTable(normMethod = "nrpm") %>%
#'   dplyr::mutate(window = paste0(seqnames, ":", start, "-", end)) %>%
#'   mesa:::getWindowAnnotation(
#'     regions = qsea::getRegions(exampleTumourNormal)[1:50],
#'     windowAnnotation = CpG_density
#'   ) %>%
#'   mesa:::makeHeatmapAnnotations(
#'     qseaSet = exampleTumourNormal,
#'     sampleAnnotation = group,
#'     windowAnnotationDf = .,
#'     useGroupMeans = FALSE
#'   ) 
#'   
makeHeatmapAnnotations <- function(qseaSet,
                                   sampleAnnotation = NULL,
                                   windowAnnotationDf = NULL,
                                   useGroupMeans = FALSE,
                                   specifiedAnnotationColors = NA,
                                   windowOrientation = "row",
                                   sampleOrientation = "column"){

  sampleAnnotationDf <- getAnnotation(qseaSet,
                                      sampleAnnotation = {{sampleAnnotation}},
                                      useGroupMeans = useGroupMeans) %>%
    dplyr::mutate_if(is.character, as.factor)

  if(is.null(windowAnnotationDf)){
    windowAnnotationDf <- data.frame()
  }

  windowAnnotationDf <- windowAnnotationDf %>%
    dplyr::mutate_if(is.character, as.factor)
  if (is.null(sampleAnnotationDf) & is.null(windowAnnotationDf)) {
    return(list(sample = NULL, window = NULL))
  }

  #Get all levels of all categorical variables and convert to color list
  levsSample <- sampleAnnotationDf %>%
    dplyr::select_if(is.factor) %>%
    purrr::map(function(x) levels(as.factor(x)))

  #Get all levels of all categorical variables and convert to color list
  levsWindow <- windowAnnotationDf %>%
    dplyr::select_if(is.factor) %>%
    purrr::map(function(x) levels(as.factor(x)))

  levs <- c(levsSample, levsWindow)

  if (length(unlist(levs)) > 0) {

    col_list_cat <- levs %>%
      unlist() %>%
      length() %>%
      hues::iwanthue(plot=FALSE, cmin = 25) %>%
      purrr::set_names(levs %>% unlist()) %>%
      utils::relist(levs) %>%
      purrr::map2(levs, purrr::set_names)
  } else {
    col_list_cat <- list()
  }

  #Make colour vectors for continuous variables
  sampleAnnotation_numeric <- sampleAnnotationDf %>%
    dplyr::select_if(is.numeric)

  #Splitting into numeric variables with all non-negative values and those with negative values
  sampleAnnotation_numeric_min_positive <- sampleAnnotation_numeric %>%
    dplyr::select_if(function(x) min(x) >= 0)

  sampleAnnotation_numeric_min_negative <- sampleAnnotation_numeric %>%
    dplyr::select_if(function(x) min(x) < 0)

  #Make colour vectors for continuous variables
  windowAnnotation_numeric <- windowAnnotationDf %>%
    dplyr::select_if(is.numeric)

  #Splitting into numeric variables with all non-negative values and those with negative values
  windowAnnotation_numeric_min_positive <- windowAnnotation_numeric %>%
    dplyr::select_if(function(x) min(x) >= 0)

  windowAnnotation_numeric_min_negative <- windowAnnotation_numeric %>%
    dplyr::select_if(function(x) min(x) < 0)

  colvecs_binary <- list(Blue = c("#fdfeff","#0053ab"),
                         Orange = c("#fffdf5","#ab5e00"),
                         Pink = c("#fffbff", "#ab0075"),
                         Green = c("#fbfff6","#36ab00"),
                         Cyan = c("#f2ffff", "#00aba2"),
                         Red = c("#fffcf9", "#ab001f"),
                         Yellow = c("#fffff4","#aba200"))

  colvecs_zerocenter <- c("BrBG","PiYG","PuOr","PRGn","RdGy") %>%
    purrr::set_names(., nm = .) %>%
    purrr::map(function(pal){
      RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal) %>%
        {c(dplyr::first(.), dplyr::nth(., 6), dplyr::last(.))}
    }
    )



  col_list_num_min_positive <- c(sampleAnnotation_numeric_min_positive, windowAnnotation_numeric_min_positive) %>%
    purrr::map2(colvecs_binary[seq_along(.)], function(val, cols){
      circlize::colorRamp2(c(min(val, na.rm = TRUE),
                             max(val, na.rm = TRUE)),
                           cols)
    }
    )

  col_list_num_min_negative <- c(sampleAnnotation_numeric_min_negative, windowAnnotation_numeric_min_negative)  %>%
    purrr::map2(colvecs_zerocenter[seq_along(.)], function(val, cols){
      circlize::colorRamp2(c(min(val, na.rm = TRUE),
                             0,
                             max(val, na.rm = TRUE)),
                           cols)
    }
    )

  annotationColors = c(col_list_cat, col_list_num_min_positive, col_list_num_min_negative)

  if(all(!is.na(specifiedAnnotationColors))){
    if(!is.list(specifiedAnnotationColors)){
      stop("Provided annotationColors object should be a list.")
    }

    commonNames <- intersect(names(annotationColors),
                             names(specifiedAnnotationColors))

    if(length(commonNames)) {
      for (i in commonNames){
        if(!all(names(annotationColors[[i]]) %in% names(specifiedAnnotationColors[[i]]))){
          missingLevels <- setdiff(names(annotationColors[[i]]), names(specifiedAnnotationColors[[i]])) %>%
            paste(collapse="\', \'")
          stop(glue::glue("Missing colors for level(s): \'{missingLevels}\' of annotation \'{i}\'."))
        }
        annotationColors[[i]] = specifiedAnnotationColors[[i]]
      }
    }

  }


  sampleAnnotation_legend_param_ls <- sampleAnnotationDf %>%
    colnames() %>%
    purrr::set_names(., nm = .) %>%
    purrr::map(function(x){
      list(name = list(direction = "horizontal"))
    })

  windowAnnotation_legend_param_ls <- windowAnnotationDf %>%
    colnames() %>%
    purrr::set_names(., nm = .) %>%
    purrr::map(function(x){
      list(name = list(direction = "horizontal"))
    })

  if(ncol(sampleAnnotationDf) > 0){
    sampleAnnot <- ComplexHeatmap::HeatmapAnnotation(which = sampleOrientation,
                                                   df    = sampleAnnotationDf,
                                                   col   = annotationColors[names(sampleAnnotationDf)],
                                                   na_col = "grey50",
                                                   annotation_legend_param = sampleAnnotation_legend_param_ls,
                                                   show_annotation_name    = FALSE)
  } else {
    sampleAnnot <- NULL
  }

  if(ncol(windowAnnotationDf) > 0){
    windowAnnot <- ComplexHeatmap::HeatmapAnnotation(which = windowOrientation,
                                                     df    = windowAnnotationDf,
                                                     col   = annotationColors[names(windowAnnotationDf)],
                                                     na_col = "grey50",
                                                     annotation_legend_param = windowAnnotation_legend_param_ls,
                                                     show_annotation_name    = FALSE)
  } else {
    windowAnnot <- NULL
  }


  return(list(sample = sampleAnnot, window = windowAnnot))
}

#'
#' Plot sample signal over windows spanning a gene (± flanks) as a heatmap using
#' **ComplexHeatmap**. The gene can be given as a HGNC symbol or Ensembl ID.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing windows and counts/betas.
#'
#' @param gene `character(1)`.  
#'   Gene identifier to plot (HGNC symbol like `"HOXA9"` or Ensembl ID like `"ENSG000..."`).
#'
#' @param normMethod `character(1)`.  
#'   Which measure to plot. One of `"beta"` or `"nrpm"`.  
#'   **Default:** `"beta"`.
#'
#' @param sampleAnnotation tidyselect specification or `character()` or `NULL`.  
#'   Columns from the sample table to display as column annotations
#'   (e.g., `c("tumour","tissue")` or bare helpers `tumour, tissue`).  
#'   **Default:** `NULL`.
#'
#' @param clusterNum `integer(1)` or `NULL`.  
#'   If set, cut the **column** dendrogram into this many clusters and display cluster labels.  
#'   **Default:** `NULL`.
#'
#' @param clusterCols `logical(1)`.  
#'   Whether to cluster columns (samples).  
#'   **Default:** `TRUE`.
#'
#' @param minDensity `numeric(1)`.  
#'   Minimum `CpG_density` required to keep a window.  
#'   **Default:** `0`.
#'
#' @param maxScale `numeric(1)`.  
#'   Upper limit of the colour scale (used for `"nrpm"`; not applied to `"beta"`).  
#'   **Default:** `1`.
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, average samples by the `group` column in the sample table (combine replicates).  
#'   **Default:** `FALSE`.
#'
#' @param upstreamDist `integer(1)`.  
#'   Number of base pairs upstream of the gene to include.  
#'   **Default:** `3000`.
#'
#' @param downstreamDist `integer(1)`.  
#'   Number of base pairs downstream of the gene to include.  
#'   **Default:** `1000`.
#'
#' @param minEnrichment `numeric(1)`.  
#'   For `"beta"`, values with enrichment `< minEnrichment` are set to `NA`.  
#'   **Default:** `3`.
#'
#' @param scaleRows `logical(1)`.  
#'   Whether to z-score scale rows (windows) before plotting.  
#'   **Default:** `FALSE`.
#'
#' @param annotationColors `list` or `NA`.  
#'   Optional named colour maps for annotations, e.g.
#'   `list(tumour = c(Tumour = "firebrick4", Normal = "blue"))`.  
#'   **Default:** `NA`.
#'
#' @param showSampleNames `logical(1)` or `NULL`.  
#'   If `NULL`, names are shown when there are fewer than 50 samples; set `TRUE`/`FALSE` to force.  
#'   **Default:** `NULL`.
#'
#' @param mart `biomaRt::Mart` or `NULL`.  
#'   Ensembl mart used to resolve gene coordinates. If `NULL`, attempts to use a mart
#'   stored on the `qseaSet` (see [setMart()]); otherwise falls back to a default for
#'   GRCh38/hg38 or hg19.  
#'   **Default:** `NULL`.
#'
#' @param idType `character(1)` or `NULL`.  
#'   Column in the mart to match `gene` against (needed for non-human/mouse setups).  
#'   **Default:** `NULL`.
#'
#' @param ... Additional arguments forwarded to **ComplexHeatmap** constructors.
#'
#' @details
#' The function resolves the gene to genomic coordinates (using `mart`, or a mart
#' stored on the object, or a built-in default for human builds), expands the
#' interval by `upstreamDist`/`downstreamDist`, filters windows by overlap and
#' `minDensity`, then draws a heatmap of `"beta"` or `"nrpm"`. For `"beta"`, values
#' below `minEnrichment` are set to `NA`. If `useGroupMeans = TRUE`, samples are
#' aggregated by `group` prior to plotting. Column annotations can be added via
#' `sampleAnnotation`; colours can be customized via `annotationColors`.
#'
#' @return
#' Draws a heatmap on the active device and (invisibly) returns the underlying
#' **ComplexHeatmap** object (e.g., a `Heatmap`/`HeatmapList`) for further composition
#' with [ComplexHeatmap::draw()].
#'
#' @seealso
#' [plotRegionsHeatmap()], [getSampleTable()], [setMart()],  
#' [biomaRt::useMart()], [ComplexHeatmap::Heatmap()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Basic gene heatmap (beta) with defaults
#' exampleTumourNormal %>% plotGeneHeatmap("HOXA9")
#'
#' # Add sample annotations (tidyselect bare names) and cluster columns into 2 groups
#' exampleTumourNormal %>%
#'   plotGeneHeatmap("HOXA9",
#'                   sampleAnnotation = c(tumour, tissue),
#'                   clusterNum = 2)
#'
#' # Custom colours and wider flanks
#' exampleTumourNormal %>%
#'   plotGeneHeatmap("HOXA9",
#'                   sampleAnnotation  = tumour,
#'                   annotationColors  = list(tumour = c(Tumour = "firebrick4", Normal = "blue")),
#'                   upstreamDist = 1000,
#'                   downstreamDist = 2000)
#'
#' 
#' @export
plotGeneHeatmap <- function(qseaSet, gene, normMethod = "beta",
                            useGroupMeans = FALSE,
                            sampleAnnotation = NULL, minDensity = 0,
                            minEnrichment = 3, maxScale = 1, clusterNum = NULL, annotationColors = NA,
                            upstreamDist = 3000, scaleRows = FALSE, clusterCols = TRUE, mart = NULL,
                            showSampleNames = NULL,
                            downstreamDist = 1000, 
                            idType = NULL,
                            ...){

  if (!is.null(getMart(qseaSet))) { mart <- getMart(qseaSet) }

  if(is.null(mart) & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"Hsapiens") & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"hg38|GRCh38")){
    mart <- biomaRt::useMart('ensembl', dataset='hsapiens_gene_ensembl', host = "https://jul2022.archive.ensembl.org")
  } else if(is.null(mart) & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"Hsapiens") & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"hg19|GRCh37")) {
    mart <- biomaRt::useMart('ensembl', dataset='hsapiens_gene_ensembl', host = "https://feb2014.archive.ensembl.org")
  } else if(is.null(mart)){
    stop("Please specify a mart object for biomaRt.")
  }

  if(is.null(idType)) {
    
    if (stringr::str_detect(gene,"^ENS.*G")) {
      idType <- "ensembl_gene_id"
    } else if (stringr::str_detect(qsea:::getGenome(qseaSet), "Hsapiens")) {
      idType <- "hgnc_symbol"
    } else if (stringr::str_detect(qsea:::getGenome(qseaSet), "Mmusculus")) {
      idType <- "mgi_symbol" 
    } else {
      stop("Please specify idType for genomes that are not human or mouse. 
           This must be a valid attribute for the given mart, see biomaRt::listAttributes.")
    }
  }


  rate <- purrr::rate_backoff(pause_base = 2, pause_min = 0.1, max_times = 3)

  # Use purrr::insistently wrapped in purrr::possibly for biomart retry logic
  safeBiomartLookup <-
    purrr::possibly(purrr::insistently(biomaRt::getBM, rate = rate, quiet = TRUE),
                    otherwise = NULL)

  bm_result <- safeBiomartLookup(
    mart = mart,
    attributes = c('hgnc_symbol', 'description', 'chromosome_name',
                   'start_position', 'end_position', 'strand', 'ensembl_gene_id'),
    filters = idType,
    values = gene
  )

  if (is.null(bm_result)) {
    stop(
      "Could not retrieve gene information from biomart after multiple attempts. Please check your internet connection or try again later."
    )
  }

  gene_details <- bm_result %>%
    dplyr::rename(seqnames = chromosome_name, start = start_position, end = end_position)

  qseaSetChr <- qseaSet %>%
    qsea::getRegions() %>%
    GenomeInfoDb::seqinfo() %>%
    GenomeInfoDb::seqnames() %>%
    stringr::str_detect("chr") %>%
    any()

  windowsChr <- gene_details %>%
    pull(seqnames) %>%
    stringr::str_detect("chr") %>%
    any()

  if (qseaSetChr & !windowsChr) {
    gene_details <- gene_details %>%
      dplyr::mutate(seqnames = paste0("chr", seqnames))
  }

  if (nrow(gene_details) > 1) {
    gene_details <- gene_details %>%
      filter(seqnames %in% GenomeInfoDb::seqlevels(qseaSet@regions))
  }

  if (nrow(gene_details) != 1) {
    stop(glue::glue("Error: {nrow(gene_details)} genes matching this name found in {mart@biomart}."))
  }

  message(glue::glue("Found {gene} on chromosome {gene_details$seqnames}, {gene_details$start} - {gene_details$end}"))
  
  gene_details_gr <- gene_details %>%
    plyranges::as_granges()

  geneGR <- gene_details %>%
    plyranges::as_granges() %>%
    plyranges::anchor_5p() %>%
    plyranges::stretch(downstreamDist) %>%
    plyranges::anchor_3p() %>%
    plyranges::stretch(upstreamDist)

  if (length(geneGR) == 0) {stop("No genomic region found!")}
  if (length(geneGR) > 1) {stop("Multiple genomic regions found!")}

  if (normMethod == "beta") {maxScale = 1}

  dataTable <- qseaSet %>%
    filterByOverlaps(regionsToOverlap = geneGR) %>%
    filterWindows(CpG_density >= minDensity) %>%
    getDataTable(normMethod = normMethod,
                 useGroupMeans = useGroupMeans,
                 minEnrichment = minEnrichment
    ) %>%
    dplyr::mutate(window = paste0(seqnames, ":",start, "-",end)) %>%
    dplyr::mutate_all( ~ dplyr::case_when(!is.nan(.x) ~ .x)) # do something with NaN values?

  if (useGroupMeans) {
    colsToFind <- qseaSet %>% getSampleGroups2() %>% names()
  } else {
    colsToFind <- qseaSet %>% qsea::getSampleNames()
  }

  numData <- dataTable %>%
    tibble::column_to_rownames("window") %>%
    dplyr::select(tidyselect::all_of(colsToFind)) %>%
    janitor::remove_empty(which = "cols", quiet = FALSE)

  if (!useGroupMeans) {

    colAnnot <- qseaSet %>%
      filter(sample_name %in% !!(colnames(numData))) %>%
      makeHeatmapAnnotations(sampleOrientation = "column",
                             useGroupMeans = useGroupMeans,
                             specifiedAnnotationColors = annotationColors,
                             sampleAnnotation = {{sampleAnnotation}} ) %>%
      purrr::pluck("sample")

  } else {
    colAnnot <- qseaSet %>%
      filter(group %in% !!(colnames(numData))) %>%
      makeHeatmapAnnotations(sampleOrientation = "column",
                             useGroupMeans = useGroupMeans,
                             specifiedAnnotationColors = annotationColors,
                             sampleAnnotation = {{sampleAnnotation}} ) %>%
      purrr::pluck("sample")
  }

  geneStrand <- gene_details_gr %>% tibble::as_tibble() %>% dplyr::pull(strand)

  annoRow <- dataTable %>%
    dplyr::select(seqnames, start, end, CpG_density, window) %>%
    plyranges::as_granges() %>%
    plyranges::mutate(annotation = dplyr::case_when(
      plyranges::count_overlaps(., gene_details_gr %>% plyranges::mutate(end = start)) > 0 & geneStrand == "+" ~ "Start",
      plyranges::count_overlaps(., gene_details_gr %>% plyranges::mutate(end = start)) > 0 & geneStrand == "-" ~ "End",
      plyranges::count_overlaps(., gene_details_gr %>% plyranges::mutate(start = end)) > 0 & geneStrand == "+" ~ "End",
      plyranges::count_overlaps(., gene_details_gr %>% plyranges::mutate(start = end)) > 0 & geneStrand == "-" ~ "Start",
      plyranges::count_overlaps(., gene_details_gr) > 0 ~ "GeneBody",
      plyranges::count_overlaps(., gene_details_gr %>% plyranges::shift_upstream(upstreamDist)) > 0 ~ "Upstream",
      plyranges::count_overlaps(., gene_details_gr %>% plyranges::shift_downstream(downstreamDist)) > 0 ~ "Downstream"
    )
    ) %>%
    as.data.frame() %>%
    dplyr::select(window, CpG_density, annotation) %>%
    tibble::column_to_rownames("window")

  windowSize <- qseaSet %>% qsea::getWindowSize()

  rowSplit <- dataTable %>%
    tibble::as_tibble() %>%
    mutate(gap = cumsum(ifelse(start - tidyr::replace_na(dplyr::lag(start),0) > !!windowSize,1,0))) %>%
    pull(gap)

  #Make row annotation object
  rowAnnot <- makeGeneHeatmapRowAnnotation(annoRow)

  #Set the cell border line width depending on the number of rows or columns, as pheatmap does.
  #Can't have gridlines for genes with too many windows, as the resulting heatmap is just grey (i.e., all you see is borders)
  if (nrow(numData) > 100 || ncol(numData) > 100) {
    rectGpParam <- grid::gpar(col = "grey", lwd = 0)
  } else {
    rectGpParam <- grid::gpar(col = "grey", lwd = 1)
  }

  #Setting a colour palette for beta-values. Could make this optional I guess
  col_fun = circlize::colorRamp2(seq(0, maxScale, length.out = 9),
                                 RColorBrewer::brewer.pal(name = "YlOrRd", n = 9))

  if (clusterCols && !is.null(clusterNum) && clusterNum > 1) {
    colSplit <- clusterNum
  } else {
    colSplit <- NULL
  }


  if (is.null(showSampleNames)) {
    if (ncol(numData) > 50) {
      showSampleNames <- FALSE
    }  else {
      showSampleNames <- TRUE
    }
  }

  annotName <- dplyr::case_when(normMethod == "beta" ~ "Beta value",
                                normMethod == "nrpm" ~ "NRPM",
                                TRUE ~ stringr::str_to_title(normMethod)
                                )

  ComplexHeatmap::Heatmap(matrix = as.matrix(numData),
                            rect_gp = rectGpParam,
                            name = annotName,
                            cluster_rows = FALSE,
                            cluster_columns = clusterCols,
                            left_annotation = rowAnnot,
                            top_annotation = colAnnot,
                            col = col_fun,
                            show_row_names = FALSE,
                            show_column_names = showSampleNames,
                            clustering_method_rows = "ward.D2",
                            clustering_method_columns = "ward.D2",
                            row_split = rowSplit,
                            row_title = NULL,
                            column_split = colSplit,
                            column_title = NULL,
                            na_col = "lightgrey",
                            heatmap_legend_param = list(legend_direction = "vertical",
                                                        at = seq(0, maxScale, length.out = 6) %>% round(1)),
                            ...
                            ) %>%
    ComplexHeatmap::draw(heatmap_legend_side = "right",
                         annotation_legend_side = "right",
                         column_title = glue::glue("{annotName} for {gene}"))

  return(invisible(numData))
}


#' Row annotations for gene heatmaps
#'
#' Helper to create a `ComplexHeatmap::HeatmapAnnotation` for **rows** showing
#' categorical tracks (e.g., region labels) and numeric tracks (e.g., CpG density).
#'
#' @param rowAnnotationDF `data.frame`/`tibble`.  
#'   One row per window with columns to annotate (e.g., `annotation`, `CpG_density`).
#'   Categorical variables (factor/character) receive discrete colour maps; numeric
#'   variables are mapped to a continuous gradient.
#'
#' @return `ComplexHeatmap::HeatmapAnnotation` configured for the row side.
#'
#' @details
#' Column types determine colour mapping automatically: factors/characters get a
#' qualitative palette; numeric columns use a continuous scale. Missing values are
#' allowed and rendered per **ComplexHeatmap** defaults.
#'
#' @seealso
#' [plotGeneHeatmap()], [plotRegionsHeatmap()], [ComplexHeatmap::HeatmapAnnotation()]
#'
#' @keywords internal
#'
#' @examples
#' set.seed(1)
#' data.frame(
#'   CpG_density = runif(10),
#'   annotation  = rep(c("Upstream","GeneBody"), 5)
#' ) %>%
#'   mesa:::makeGeneHeatmapRowAnnotation() 
#'   
makeGeneHeatmapRowAnnotation <- function(rowAnnotationDF){

  annotationColDf <- rowAnnotationDF %>%
    dplyr::mutate_if(is.character, as.factor)

  #Get all levels of all categorical variables and convert to color list
  levs <- annotationColDf %>%
    dplyr::select_if(is.factor) %>%
    purrr::map(function(x) levels(as.factor(x)))

  col_list_cat <- levs %>%
    unlist() %>%
    length() %>%
    hues::iwanthue(plot=FALSE) %>%
    purrr::set_names(levs %>% unlist()) %>%
    utils::relist(levs) %>%
    purrr::map2(levs, purrr::set_names)

  #Make colour vectors for continuous variables
  annotationCol_numeric <- annotationColDf %>%
    dplyr::select_if(is.numeric)

  #Splitting into numeric variables with all non-negative values and those with negative values
  annotationCol_numeric_min_positive <- annotationCol_numeric %>%
    dplyr::select_if(function(x) min(x) >= 0)

  annotationCol_numeric_min_negative <- annotationCol_numeric %>%
    dplyr::select_if(function(x) min(x) < 0)

  colvecs_binary <- c("Reds","YlGnBu","YlOrBr","PuRd","Blues","Purples") %>%
      purrr::set_names(., nm = .) %>%
      purrr::map(function(pal){
        RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal) %>%
          {c(dplyr::first(.), dplyr::last(.))}
        }
    )

  colvecs_zerocenter <- c("BrBG","PiYG","PuOr","PRGn","RdGy") %>%
    purrr::set_names(., nm = .) %>%
    purrr::map(function(pal){
      RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal) %>%
        {c(dplyr::first(.), dplyr::nth(., 6), dplyr::last(.))}
      }
    )

  col_list_num_min_positive <- annotationCol_numeric_min_positive %>%
    purrr::map2(colvecs_binary[seq_along(ncol(.))], function(val, cols){
      circlize::colorRamp2(c(min(val, na.rm = TRUE),
                             max(val, na.rm = TRUE)),
                           cols)
      }
    )

  col_list_num_min_negative <- annotationCol_numeric_min_negative %>%
    purrr::map2(colvecs_zerocenter[seq_along(ncol(.))], function(val, cols){
      circlize::colorRamp2(c(min(val, na.rm = TRUE),
                             0,
                             max(val, na.rm = TRUE)),
                           cols)
      }
    )

  annotationColors = c(col_list_cat, col_list_num_min_positive, col_list_num_min_negative)

  annotation_legend_param_ls <- annotationColDf %>%
    colnames() %>%
    purrr::set_names(., nm = .) %>%
    purrr::map(function(x){
      list(name = list(direction = "horizontal"))
      }
    )

  annot <- ComplexHeatmap::HeatmapAnnotation(which = "row",
                                             df    = annotationColDf,
                                             col   = annotationColors,
                                             annotation_legend_param = annotation_legend_param_ls,
                                             show_annotation_name    = FALSE)

  return(annot)
}


#' Distribution of windows across genomic features
#'
#' Annotate windows (e.g., promoters/exons/introns) and plot the distribution
#' per sample as stacked/dodged/filled bars.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing windows and signal (beta or nrpm).
#'
#' @param cutoff `numeric(1)`.  
#'   Threshold applied to the chosen normalisation measure per window.  
#'   **Default:** `1`.
#'
#' @param barType `character(1)`.  
#'   Bar layout: one of `"stack"`, `"dodge"`, or `"fill"` (relative proportions).  
#'   **Default:** `"stack"`.
#'
#' @param normMethod `character(1)`.  
#'   Normalisation/measure to threshold. One of `"nrpm"` or `"beta"`.  
#'   **Default:** `"nrpm"`.
#'
#' @return A `ggplot2` object showing counts (or proportions, for `"fill"`)
#'   of feature classes per sample above `cutoff`.
#'
#' @details
#' Feature classes are taken from region annotations available on the windows
#' (e.g., columns produced by [annotateWindows()] such as `shortAnno`/`annotation`,
#' if present). Bars are positioned according to `barType` (`stack`/`dodge`/`fill`).
#'
#' @seealso
#' [annotateWindows()], [qsea::makeTable()], [ggplot2::geom_bar()]
#' 
#' @family annotation-summaries
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Beta >= 0.75, stacked bars per sample
#' exampleTumourNormal %>%
#'   plotGenomicFeatureDistribution(normMethod = "beta", cutoff = 0.75)
#'
#' # NRPM >= 2, show within-sample proportions
#' exampleTumourNormal %>%
#'   plotGenomicFeatureDistribution(normMethod = "nrpm", cutoff = 2, barType = "fill")
#'   
#' @export
plotGenomicFeatureDistribution <- function(qseaSet, cutoff = 1 , barType = "stack", normMethod = "nrpm"){

  #TODO: Rewrite this function to work with any genome! Needs more options exposed (TxDb etc).
  
  temp <- qseaSet %>%
    qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = normMethod) %>%
    # makeTable(keep = which(qsea::getCounts(.)[,sampleName] >= (0.75 * rpmFactor)), samples = sampleName, norm_methods = "nrpm") %>%
    # filter(!!dplyr::sym(paste0(sampleName,"_nrpm")) >= cutoff) %>%
    qseaTableToChrGRanges() %>%
    ChIPseeker::annotatePeak(tssRegion = c(-2000, 500),
                             level = "gene",
                             TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                             annoDb = "org.Hs.eg.db",
                             overlap = "all",
                             verbose = FALSE)

  # temp2 <- temp@anno %>%
  #   group_by(annoShort)

  featureTable <- purrr::map_dfr(qsea::getSampleNames(qseaSet),function(x){
    temp@anno %>%
      dplyr::filter(!!dplyr::sym(paste0(x,"_",normMethod)) > cutoff) %>%
      dplyr::mutate(annoShort = stringr::str_replace(annotation, "on \\(.*", "on")) %>%
      tibble::as_tibble() %>%
      dplyr::pull(annoShort) %>%
      table() %>%
      tibble::enframe(name = "feature") %>%
      dplyr::mutate(sample = x)
  }
  )

  featureTable %>%
    ggplot2::ggplot(ggplot2::aes(y = value, x = sample, fill = feature)) +
    ggplot2::geom_bar(position = barType, stat = "identity") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    ggplot2::labs(x = "Sample",
                  y = "Fraction",
                  legend = "Feature",
                  subtitle = glue::glue("{getWindowSize(qseaSet)}bp windows with at least {cutoff} {normMethod}"))

}

#' Sample correlation heatmap
#'
#' Compute and plot the correlation matrix across samples (or group means)
#' using the selected normalisation measure and optional region filters.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing windows and signal (beta or nrpm).
#'
#' @param regionsToOverlap `GRanges`, `data.frame` (coercible to `GRanges`), or `NULL`.  
#'   If provided, restrict computation to windows overlapping these regions; otherwise
#'   use all windows.  
#'   **Default:** `NULL`.
#'
#' @param useGroupMeans `logical(1)`.  
#'   If `TRUE`, average replicates by the `group` column and compute correlations
#'   on group means; if `FALSE`, use individual samples.  
#'   **Default:** `FALSE`.
#'
#' @param sampleAnnotation tidyselect specification or `character()` or `NULL`.  
#'   Columns from the sample table to display as annotations on the heatmap
#'   (e.g., `c("tumour","patient")` or bare helpers `tumour, patient`).  
#'   **Default:** `NULL`.
#'
#' @param normMethod `character(1)`.  
#'   Measure to correlate. One of `"nrpm"` or `"beta"`.  
#'   **Default:** `"nrpm"`.
#'
#' @param minDensity `numeric(1)`.  
#'   Minimum `CpG_density` required to keep a window.  
#'   **Default:** `0`.
#'
#' @param minEnrichment `numeric(1)`.  
#'   For `"beta"`, values with enrichment `< minEnrichment` are set to `NA` before
#'   correlation.  
#'   **Default:** `3`.
#'
#' @param annotationColors `list` or `NA`.  
#'   Optional named colour maps for annotations (passed to **pheatmap**), e.g.
#'   `list(tumour = c(Tumour = "firebrick4", Normal = "blue"))`.  
#'   **Default:** `NA`.
#'
#' @param ...  
#'   Additional arguments forwarded to [pheatmap::pheatmap()].
#'
#' @details
#' Windows are optionally restricted by `regionsToOverlap` and filtered by `minDensity`.
#' Correlations are computed on the chosen `normMethod` after any group averaging.
#' For `"beta"`, entries below `minEnrichment` are set to `NA` to down-weight
#' low-enrichment windows.
#'
#' @return A [`pheatmap::pheatmap`] object.
#'
#' @seealso
#' [plotRegionsHeatmap()], [getDataTable()], [pheatmap::pheatmap()]
#'
#' @family heatmaps
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # By default uses NRPM
#' exampleTumourNormal %>% plotCorrelationMatrix()
#'
#' # Using beta values and adding sample annotations
#' exampleTumourNormal %>% 
#'   plotCorrelationMatrix(normMethod = "beta",
#'                         sampleAnnotation = c(tumour, patient))
#'
#' @export
plotCorrelationMatrix <- function(qseaSet, regionsToOverlap = NULL, useGroupMeans = FALSE, sampleAnnotation = NULL, normMethod = "nrpm",
                                  minEnrichment = 3, annotationColors = NA, minDensity = 0, ...){

  ##TODO: Swap from pheatmap to ComplexHeatmap 
  if (!is.null(regionsToOverlap)) {
    regionsToOverlap <- asValidGranges(regionsToOverlap)

    if (length(regionsToOverlap) == 0) {stop("No genomic regions given!")}

    qseaSet <- qseaSet %>%
      filterByOverlaps(regionsToOverlap = regionsToOverlap)
  }

  annotationDf = getAnnotation(qseaSet, sampleAnnotation = {{sampleAnnotation}}, useGroupMeans = useGroupMeans)

  if (ncol(annotationDf) == 0) {
    annotationDf <- NULL
  }
  
  dataTab <- qseaSet %>%
    filterWindows(CpG_density >= minDensity) %>%
    getDataTable(normMethod = normMethod,
                 useGroupMeans = useGroupMeans,
                 minEnrichment = minEnrichment
    ) %>%
    dplyr::mutate(window = paste0(seqnames, ":",start, "-",end)) %>%
    dplyr::mutate_all( ~ dplyr::case_when(!is.nan(.x) ~ .x)) # do something with NaN values?

  if(useGroupMeans){
    colsToFind <- qseaSet %>% getSampleGroups2() %>% names()
  } else {
    colsToFind <- qseaSet %>% qsea::getSampleNames()
  }

  numData <- dataTab %>%
    dplyr::select(tidyselect::all_of(colsToFind)) %>%
    janitor::remove_empty(which = "cols", quiet = FALSE) %>%
    janitor::remove_empty(which = "rows", quiet = FALSE)

  corMat <- numData %>%
    stats::cor(use = "pairwise.complete.obs") %>%
    as.data.frame()
  
  corMat %>%
    pheatmap::pheatmap(display_numbers = TRUE,
                       color = RColorBrewer::brewer.pal(name = "YlOrRd", n = 9),
                       annotation_row = annotationDf,
                       annotation_col = annotationDf,
                       annotation_colors = annotationColors, ...)
  
  return(invisible(corMat))
}


#' UpSet plot of DMR overlaps
#'
#' Visualise the overlap of significant DMR sets across contrasts (columns ending
#' with `*_adjPval`) using an UpSet plot.
#'
#' @param DMRtable `data.frame`.  
#'   Table as returned by [calculateDMRs()] (optionally pre-filtered).
#'
#' @param string `character(1)` or `NULL`.  
#'   Optional regular expression to subset the set/contrast names (applied after
#'   stripping suffixes).  
#'   **Default:** `NULL` (use all sets).
#'
#' @param removeVS `logical(1)`.  
#'   If `TRUE`, remove the `"_vs_"` substring and everything after it from set
#'   names (e.g., `"Tumour_vs_Normal"` → `"Tumour"`).  
#'   **Default:** `FALSE`.
#'
#' @param minAdjPval `numeric(1)`.  
#'   Adjusted P-value threshold; windows with `adjPval <= minAdjPval` are included
#'   in each set.  
#'   **Default:** `0.05`.
#'
#' @param ...  
#'   Additional arguments forwarded to [UpSetR::upset()].
#'
#' @details
#' The function discovers DMR sets by locating columns that end with `"_adjPval"`.
#' For each such column, a logical membership is formed at the chosen FDR cutoff
#' (`minAdjPval`). Optionally, set names are simplified with `removeVS`, and then
#' subset with `string` if provided. The resulting membership matrix is visualised
#' with **UpSetR**.
#'
#' @return
#' Draws an UpSet plot on the active device and **invisibly** returns the result
#' from [UpSetR::upset()] (for further customization if needed).
#'
#' @seealso
#' [calculateDMRs()], [UpSetR::upset()]
#'
#' @family dmr-plots
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Overlap of significant DMR sets across contrasts (default FDR 0.05)
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type", contrasts = "all_vs_NormalLung") %>%
#'   plotDMRUpset()
#'
#' # Clean set names by dropping the '_vs_*' suffix
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type", contrasts = "all_vs_NormalLung") %>%
#'   plotDMRUpset(removeVS = TRUE)
#'
#' @export
plotDMRUpset <- function(DMRtable, string = NULL, removeVS = FALSE, minAdjPval = 0.05, ...){

  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop(
      "Package \"UpSetR\" must be installed to use this function.",
      call. = FALSE
    )
  }

  temp <- DMRtable %>%
    dplyr::select(tidyselect::matches("adjPval")) %>%
    {if (removeVS) {dplyr::rename_with(., ~stringr::str_remove(.,"_vs_.*")) } else {dplyr::rename_with(., ~stringr::str_remove(.,"_adjPval$"))}} %>%
    {if (!is.null(string)) {dplyr::select(., tidyselect::matches(string))} else . }

  if (ncol(temp) < 2) {stop(glue::glue("Only {ncol(temp)} column remaining, nothing to plot!"))}

  purrr::map(stats::setNames(colnames(temp),colnames(temp)),function(x){  which(temp[,x,drop = FALSE] <= minAdjPval)}) %>%
    UpSetR::fromList() %>%
    UpSetR::upset(nsets = ncol(.), nintersects = 100, order.by = "freq", text.scale	 = 1.8, ...)
}


#' Sample/group annotation data.frame
#'
#' Note that sampleAnnotation must be enclosed within double curly brackets when used.
#'
#' Select sample-table columns to use as annotations, either per-sample or
#' per-group (with consistency checks across replicates).
#'
#' @param qseaSet `qseaSet`.  
#'   Object from which sample metadata are taken.
#'
#' @param useGroupMeans `logical(1)`.  
#'   Build annotations at the **group** level (aggregate replicates by `group`);
#'   if `FALSE`, return per-sample annotations.  
#'   **Default:** `FALSE`.
#'
#' @param sampleAnnotation tidyselect specification, `character()`, or `NULL`.  
#'   Columns from the sample table to include as annotations. Accepts unquoted
#'   tidyselect (e.g., `tumour, type`) or a character vector (e.g., `c("tumour","type")`).  
#'   **Default:** `NULL` (no extra columns beyond the identifier).
#'
#' @return `data.frame` of annotations with row names set to `sample_name`
#'   (or to `group` when `useGroupMeans = TRUE`).
#'
#' @details
#' When *programming* with tidyselect inside your own functions, wrap the
#' `sampleAnnotation` argument in **double curly braces** (`{{ }}`) to enable
#' tidy-evaluation. End-users can pass unquoted column names or character
#' vectors directly.
#'
#' @seealso
#' [makeHeatmapAnnotations()], [plotRegionsHeatmap()], [plotCorrelationMatrix()]
#'
#' @family heatmap-annotation
#'
#' @keywords internal
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Works with none, one, or multiple columns (unquoted tidyselect):
#' exampleTumourNormal %>% mesa:::getAnnotation()
#' exampleTumourNormal %>% mesa:::getAnnotation(sampleAnnotation = type)
#' exampleTumourNormal %>% mesa:::getAnnotation(sampleAnnotation = c(tumour, type))
#'
#' # Also works with quoted column names:
#' exampleTumourNormal %>% mesa:::getAnnotation(sampleAnnotation = "type")
#' exampleTumourNormal %>% mesa:::getAnnotation(sampleAnnotation = c("tumour","type"))
#'
#' # group-level annotations for tumour & tissue
#' exampleTumourNormal %>%
#'   dplyr::mutate(group = stringr::str_remove(sample_name, "[0-9]")) %>%
#'   mesa:::getAnnotation(sampleAnnotation = c(tumour, tissue), useGroupMeans = TRUE)
#'
#' exampleTumourNormal %>%
#'   dplyr::mutate(group = stringr::str_remove(sample_name, "[0-9]")) %>%
#'   mesa:::getAnnotation(sampleAnnotation = tumour, useGroupMeans = TRUE)
#'   
getAnnotation <- function(qseaSet, useGroupMeans = FALSE, sampleAnnotation = NULL){

  if (rlang::quo_is_null(rlang::enquo(sampleAnnotation))) {
    return(data.frame())
  }

  annotationColDf <- data.frame()

  if (!useGroupMeans) {
    annotationColDf <- qseaSet %>%
      qsea::getSampleTable() %>%
      #dplyr::arrange(sample_name) %>%
      dplyr::select(!!!rlang::enquos(sampleAnnotation))

  } else if (useGroupMeans) {
    groupSampleTab <- qseaSet %>%
      qsea::getSampleTable() %>%
      dplyr::select(group, !!!rlang::enquos(sampleAnnotation))

    distinctGroupTab <- groupSampleTab %>%
      tibble::remove_rownames() %>%
      dplyr::distinct()

    if (nrow(distinctGroupTab) != (qseaSet %>% getSampleGroups2() %>% length())) {
      problemGroups <- distinctGroupTab %>%
        dplyr::count(group) %>%
        filter(n > 1) %>%
        pull(group) %>%
        paste(., collapse = "; ")
      stop(glue::glue("Grouped annotation contains differing annotation in {problemGroups}"))
    }

    annotationColDf <- distinctGroupTab %>%
      #arrange(group) %>%
      tibble::column_to_rownames("group") %>%
      as.data.frame()
  }

  if (ncol(annotationColDf) == 0) {
    annotationColDf <- data.frame()
  }

  return(annotationColDf)
}


#' Test wrapper for \code{plotGeneHeatmap} that skips on network errors
#'
#' Calls \code{plotGeneHeatmap(...)} inside \code{tryCatch}. If a likely
#' BiomaRt/network error occurs, the test is skipped instead of failing.
#'
#' @param ... Arguments passed to \code{plotGeneHeatmap()}.
#'
#' @details Errors are matched case-insensitively against patterns such as
#'   "biomart", "SSL", "connection", "timeout", "could not resolve",
#'   "unexpected eof", and "http 500". On a match, the test is skipped via
#'   \code{testthat::skip()}.
#'
#' @return Invisibly returns \code{TRUE} when the plot call succeeds and the
#'   test is marked as succeeded. Otherwise, rethrows non-network errors.
#'
#' @seealso \code{\link{plotGeneHeatmap}}
#'
#' @keywords internal
testPlotGeneHeatmap <- function(...) {
  tryCatch({
    plotGeneHeatmap(...)
    testthat::succeed()
  }, error = function(e) {
    if (stringr::str_detect(e$message, "(?i)biomart|SSL|connection|timeout|could not resolve|unexpected eof|http 500")) {
      testthat::skip("Connection to biomart failed, skipping test.")
    } else {
      stop(e)
    }
  })
}