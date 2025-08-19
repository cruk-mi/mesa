#' Heatmap of signal across genomic regions
#'
#' Plot a heatmap of per-window values (NRPM/beta) for regions of interest,
#' optionally adding sample and window annotations and clustering.
#'
#' @param qseaSet qsea::qseaSet.
#' @param regionsToOverlap GenomicRanges::GRanges or data frame coercible to GRanges;
#'   if `NULL`, uses all regions in `qseaSet`.
#' @param normMethod character(1) normalisation method to plot, typically
#'   `"nrpm"` or `"beta"`. Defaults to `"beta"`.
#' @param useGroupMeans logical(1) whether to average replicates by the `group`
#'   column of the sample table (i.e. combine replicates)
#' @param sampleAnnotation tidyselect-style columns from the sample table to show
#'   as annotations (use unquoted names, e.g. `group`).
#' @param windowAnnotation tidyselect-style columns from `qseaSet` regions or
#'   `regionsToOverlap` to show as row annotations (unquoted).
#' @param clusterNum integer(1) optional number of column clusters (if
#'   `clusterCols = TRUE`) to break the column dendrogram into.
#' @param maxScale numeric(1) upper bound of colour scale (ignored for `"beta"`).
#' @param clip numeric(1) clip values above this threshold before plotting.
#' @param minDensity numeric(1) minimum CpG density for windows to be kept.
#' @param clusterRows,clusterCols logical whether to cluster rows/columns.
#' @param clusterMethod character(1) clustering method passed to ComplexHeatmap.
#' @param annotationColors optional `list` of colour mappings for annotation
#'   variables; if omitted, colours are generated.
#' @param minEnrichment integer(1) minimum reads for non-NA beta (qsea rule).
#' @param showSampleNames logical(1) whether to draw sample names (default
#'   `TRUE` if ≤ 50 samples).
#' @param annotationPosition character(1) legend side (`"right"`, `"bottom"`).
#' @param title optional character(1) column title.
#' @param ... additional arguments passed to `ComplexHeatmap::Heatmap()`.
#' 
#' @return Invisibly returns the numeric matrix used for plotting; the heatmap
#'   is drawn to the active graphics device.
#' @seealso [getDataTable()], [getWindowAnnotation()], [makeHeatmapAnnotations()],
#'   [plotGeneHeatmap()]
#' @family heatmaps
#'
#' @examples
#' \donttest{
#' # Plot beta values for a small subset of windows with sample annotations
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' regs <- qsea::getRegions(qs)[1:400]   # small region subset for speed
#' plotRegionsHeatmap(qs, regs, normMethod = "beta",
#'                    sampleAnnotation = group,              # add sample group
#'                    windowAnnotation = CpG_density,        # row-side CpG info
#'                    clusterCols = TRUE, clusterRows = FALSE)
#' }
#' 
#' @export
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
                            left_annotation = rowAnnot) %>%
    ComplexHeatmap::draw(heatmap_legend_side = "bottom",
                         annotation_legend_side = annotationPosition,
                         column_title = title,
                         column_title_gp = grid::gpar(fontsize = 16))

}


#' Build row annotations for region heatmaps
#'
#' Prepare a data.frame of window-level annotations aligned to rows in the
#' heatmap created by [plotRegionsHeatmap()].
#'
#' @param dataTab tibble/data.frame produced by [getDataTable()] plus `seqnames`,
#'   `start`, `end`, and a `window` column (`chr:start-end`).
#' @param regions GenomicRanges::GRanges for overlap-derived annotations
#'   (e.g., columns already present on `regions`).
#' @param windowAnnotation tidyselect-style columns (unquoted) to keep as
#'   annotation (from either `dataTab` or `regions`).
#' @param clusterRows logical; if `FALSE`, preserves genomic order of windows.
#'
#' @return `data.frame` with row names set to the `window` identifiers and the
#'   requested annotation columns.
#' @seealso [plotRegionsHeatmap()], [makeHeatmapAnnotations()]
#' @keywords internal
#'
#' @examples
#' # Build a minimal window annotation table from a toy selection
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' regs <- qsea::getRegions(qs)[1:50]
#' dt   <- qs |>
#'   filterByOverlaps(regs) |>
#'   getDataTable(normMethod = "nrpm") |>
#'   dplyr::mutate(window = paste0(seqnames, ":", start, "-", end))
#' wa <- getWindowAnnotation(dt, regs, windowAnnotation = CpG_density)
#' head(wa)
#' }
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

  rowAnnotDf %>%
    dplyr::select(-.rowID) %>%
    tibble::column_to_rownames("window")

}


#' Create ComplexHeatmap annotation objects (samples and windows)
#'
#' Helper to generate `HeatmapAnnotation`s for samples and windows with
#' automatic colour handling for categorical and numeric variables.
#'
#' @param qseaSet qsea::qseaSet.
#' @param sampleAnnotation tidyselect-style columns from the sample table
#'   to annotate (unquoted).
#' @param windowAnnotationDf data.frame of row annotations as returned by
#'   [getWindowAnnotation()].
#' @param useGroupMeans logical whether to annotate groups instead of samples.
#' @param specifiedAnnotationColors optional `list` mapping levels → colours to
#'   override the auto-generated palette.
#' @param windowOrientation character(1) orientation for window annotation
#'   (`"row"` or `"column"`).
#' @param sampleOrientation character(1) orientation for sample annotation
#'   (`"column"` or `"row"`).
#'
#' @return `list(sample = <HeatmapAnnotation|NULL>, window = <HeatmapAnnotation|NULL>)`
#'   ready to pass to `ComplexHeatmap::Heatmap(...)`.
#' @seealso [plotRegionsHeatmap()], [getWindowAnnotation()]
#' @keywords internal
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' wa <- data.frame()  # no window annotations
#' anns <- makeHeatmapAnnotations(qs, sampleAnnotation = group,
#'                                windowAnnotationDf = wa, useGroupMeans = FALSE)
#' str(anns)
#' }
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


#' Heatmap across a single gene locus
#'
#' Retrieve a gene locus via biomaRt and plot window-level values (NRPM/beta)
#' across upstream/gene body/downstream regions with sample annotations.
#'
#' @param qseaSet qsea::qseaSet.
#' @param gene character(1) HGNC/MGI symbol or Ensembl gene ID.
#' @param normMethod character(1) `"nrpm"` or `"beta"`.
#' @param useGroupMeans logical collapse replicates by `group`.
#' @param sampleAnnotation tidyselect-style sample-table columns to annotate.
#' @param minDensity numeric minimum CpG density.
#' @param minEnrichment integer minimum reads for non-NA beta.
#' @param maxScale numeric upper colour limit (ignored for `"beta"`).
#' @param clusterNum integer optional number of column clusters.
#' @param annotationColors optional list of colours for annotations.
#' @param upstreamDist,downstreamDist integer(1) bp to extend around gene.
#' @param scaleRows logical whether to z-scale rows (not used here; reserved).
#' @param clusterCols logical cluster columns.
#' @param mart optional biomaRt Mart; if `NULL`, a reasonable default is chosen
#'   for human (hg38/hg19) otherwise `idType` must be provided.
#' @param showSampleNames logical show sample names (default `TRUE` if ≤ 50).
#' @param idType character(1) biomaRt attribute for `gene` identifiers (e.g.
#'   `"ensembl_gene_id"`, `"hgnc_symbol"`, `"mgi_symbol"`).
#' @param ... passed to `ComplexHeatmap::Heatmap`.
#'
#' @return Invisibly returns the numeric matrix used for plotting; the heatmap
#'   is drawn to the graphics device.
#' @seealso [plotRegionsHeatmap()], [makeHeatmapAnnotations()]
#' @family heatmaps
#'
#' @examples
#' \donttest{
#' # Plot around a gene symbol (human dataset)
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' plotGeneHeatmap(qs, gene = "TP53", normMethod = "beta",
#'                 upstreamDist = 3000, downstreamDist = 1000,
#'                 sampleAnnotation = group, clusterCols = TRUE)
#' }
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

  gene_details <- biomaRt::getBM(mart = mart,
                                 attributes = c('hgnc_symbol', 'description', 'chromosome_name',
                                                'start_position', 'end_position', 'strand','ensembl_gene_id'),
                                 filters = idType,
                                 values = gene) %>%
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

  #Make rowannoation object
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
                                                        at = seq(0, maxScale, length.out = 6) %>% round(1))
                            ) %>%
    ComplexHeatmap::draw(heatmap_legend_side = "right",
                         annotation_legend_side = "right",
                         column_title = glue::glue("{annotName} for {gene}"))

  return(invisible(numData))
}


#' Row annotations for gene heatmaps
#'
#' Helper to create a `HeatmapAnnotation` for rows showing categorical
#' tracks (e.g., region labels) and numeric tracks (e.g., CpG density).
#'
#' @param rowAnnotationDF data.frame with row-wise annotation variables
#'   (e.g., columns `annotation`, `CpG_density`) indexed by window.
#'
#' @return `ComplexHeatmap::HeatmapAnnotation` for the row side.
#' @keywords internal
#'
#' @examples
#' \donttest{
#' df <- data.frame(CpG_density = runif(10), annotation = rep(c("Upstream","GeneBody"),5))
#' makeGeneHeatmapRowAnnotation(df)
#' }
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

  #if(is.null())

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
    purrr::map2(colvecs_binary[1:ncol(.)], function(val, cols){
      circlize::colorRamp2(c(min(val, na.rm = TRUE),
                             max(val, na.rm = TRUE)),
                           cols)
      }
    )

  col_list_num_min_negative <- annotationCol_numeric_min_negative %>%
    purrr::map2(colvecs_zerocenter[1:ncol(.)], function(val, cols){
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
#' @param qseaSet qsea::qseaSet.
#' @param cutoff numeric threshold on the chosen `normMethod`.
#' @param barType character(1) `"stack"`, `"dodge"`, or `"fill"`.
#' @param normMethod character(1) `"nrpm"` or `"beta"`.
#'
#' @return A `ggplot` object of the distribution
#' @seealso [getGenomicFeatureDistribution()] (tabular summary), [annotateWindows()]
#' @family annotation-summaries
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' plotGenomicFeatureDistribution(qs, cutoff = 1, barType = "fill", normMethod = "nrpm")
#' }
#' 
#' @export
plotGenomicFeatureDistribution <- function(qseaSet, cutoff = 1 , barType = "stack", normMethod = "nrpm"){

  #rpmFactor <- (qsea::getLibSize(qseaSetcfDNA)[sampleName]/1000000)

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
#' using the selected normalisation method and optional region filters.
#'
#' @param qseaSet qsea::qseaSet.
#' @param regionsToOverlap optional GRanges to restrict windows.
#' @param useGroupMeans logical average replicates by `group`.
#' @param sampleAnnotation tidyselect-style sample-table columns to annotate.
#' @param normMethod character(1) `"nrpm"` or `"beta"`.
#' @param minDensity numeric minimum CpG density to keep windows.
#' @param minEnrichment integer minimum reads for non-NA beta.
#' @param annotationColors optional list of colours for annotations (pheatmap).
#' @param ... passed to `pheatmap::pheatmap()`.
#'
#' @return A `pheatmap` object.
#' @seealso [plotRegionsHeatmap()], [getDataTable()]
#' @family heatmaps
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' plotCorrelationMatrix(qs, normMethod = "nrpm", sampleAnnotation = group)
#' }
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

  numData %>%
    stats::cor(use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    pheatmap::pheatmap(display_numbers = TRUE,
                       color = RColorBrewer::brewer.pal(name = "YlOrRd", n = 9),
                       annotation_row = annotationDf,
                       annotation_col = annotationDf,
                       annotation_colors = annotationColors, ...)
}


#' UpSet plot of DMR overlaps
#'
#' Visualise overlap of significant DMR sets across contrasts from a DMR table
#' (columns `*_adjPval`) using an UpSet plot.
#'
#' @param DMRtable data.frame as returned by [calculateDMRs()] (possibly filtered).
#' @param string optional regex to subset DMR columns (after stripping suffixes).
#' @param removeVS logical remove `"_vs_"` and following text from column names.
#' @param minAdjPval numeric adjusted P-value threshold for inclusion.
#' @param ... passed to `UpSetR::upset()`.
#'
#' @return An UpSet plot (drawn), and the function returns the `UpSetR` object.
#' @seealso [calculateDMRs()]
#' @family dmr-plots
#'
#' @examples
#' \donttest{
#' # Suppose `dmr` is a DMR table from calculateDMRs()
#' # plotDMRUpset(dmr, string = "LUAD|LUSC", minAdjPval = 0.05)
#' }
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
#' Select sample-table columns to use as annotations, either per-sample or
#' per-group (with consistency checks across replicates). Note that sampleAnnotation 
#' must be enclosed within double curly brackets when used.
#'
#' @param qseaSet qsea::qseaSet.
#' @param useGroupMeans logical build annotations at group level.
#' @param sampleAnnotation tidyselect-style columns to include (unquoted).
#'
#' @return data.frame of annotations; rows named by `sample_name` (or `group`).
#' @seealso [makeHeatmapAnnotations()], [plotRegionsHeatmap()], [plotCorrelationMatrix()]
#' @family heatmap-annotation
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- exampleTumourNormal
#' ann <- getAnnotation(qs, useGroupMeans = FALSE, sampleAnnotation = group)
#' head(ann)
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
