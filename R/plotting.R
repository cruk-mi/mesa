
#' This function takes a qseaSet and a GRanges object, and plots the expression across the gene as a heatmap
#' @param qseaSet The qseaSet object.
#' @param regionsToOverlap A genomic ranges to plot.
#' @param normMethod Whether to plot nrpm values or beta values.
#' @param useGroupMeans Whether to average samples over the group column (i.e. combine replicates)
#' @param sampleAnnotation A character vector with names of columns from the sampleTable to use as annotation, or a pre-defined annotation data frame.
#' @param clusterNum A number of clusters to break the column dendrogram into.
#' @param maxScale The maximum of the scale, not used when plotting beta values.
#' @param clip An upper level to clip the data at, anything higher than this is replaced with the clip value.
#' @param minDensity A minimum CpG density level to filter out windows with values lower than.
#' @param clusterRows Whether to cluster the rows or not.
#' @param clusterCols Whether to cluster the columns or not.
#' @param clusterMethod What method to use to cluster the dendrograms.
#' @param annotationColors A list specifying some or all of the colours to use for the annotations.
#' @param minEnrichment Minimum enrichment factor for beta values, will give NAs below this.
#' @param annotationPosition Where to put the annotation for the samples, e.g. "bottom" or "right".
#' @param title A title to add to the top of the plot.
#' @param ... Other arguments to pass to ComplexHeatmap.
#' @return A qseaSet object with the sampleTable enhanced with the information on number of reads etc
#' @export
plotRegionsHeatmap <- function(qseaSet, regionsToOverlap = NULL,
                                normMethod = "beta",
                                sampleAnnotation = NULL,
                                annotationColors = NA,
                                useGroupMeans = FALSE,
                                clusterRows = FALSE,
                                clusterCols = TRUE,
                                minEnrichment = 3,
                                maxScale = 5,
                                clusterNum = 2,
                                clip = 1000000000,
                                minDensity = 0,
                                annotationPosition = "right",
                                title = NULL,
                                clusterMethod = "ward.D2", ...) {

  #build the column annotation. Need the {{ }} to parse either strings or tidy selections
  colAnnot <- makeHeatmapAnnotation(qseaSet,
                                    orientation = "column",
                                    useGroupMeans = useGroupMeans,
                                    sampleAnnotation = {{sampleAnnotation}} )

  if (is.null(regionsToOverlap)) {
    regionsToOverlap <- qseaSet %>%
      qsea::getRegions()
  }

  regionsToOverlap <- asValidGranges(regionsToOverlap)

  if (length(regionsToOverlap) > 20000) {stop("More than 20000 regions requested.")}

  if (length(regionsToOverlap) == 0) {stop("No genomic regions given!")}

  if (normMethod == "beta") {maxScale = min(clip,1)}

  colour_palette <- RColorBrewer::brewer.pal(name = "YlOrRd", n = 9)

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
    colsToFind <- qseaSet %>% qsea::getSampleGroups() %>% names()
  } else {
    colsToFind <- qseaSet %>% qsea::getSampleNames()
  }

  numData <- dataTab %>%
    dplyr::select(tidyselect::all_of(colsToFind)) %>%
    clipFn(a = 0, b = clip) %>%
    janitor::remove_empty(which = "cols", quiet = FALSE) %>%
    janitor::remove_empty(which = "rows", quiet = FALSE)

  if (ncol(dataTab) == 1) {
    clusterCols <- FALSE
  }

  if (clusterRows) {
    dataTab <- remove_almost_empty_rows(dataTab)
  }

  if (clusterCols) {
    colSplit <- clusterNum
  } else {
    colSplit <- NULL
  }


  numData %>%
    as.matrix() %>%
    ComplexHeatmap::Heatmap(name = "Beta value",
                            cluster_rows = clusterRows,
                            cluster_columns = clusterCols,
                            show_row_names = FALSE,
                            col = colour_palette,
                            clustering_method_rows = clusterMethod,
                            clustering_method_columns = clusterMethod,
                            heatmap_legend_param = list(legend_direction = "horizontal",
                                                        at = seq(0, maxScale, length.out = 6) %>% round(1)),
                            column_split = colSplit,
                            column_title = NULL,
                            top_annotation = colAnnot) %>%
    ComplexHeatmap::draw(heatmap_legend_side = "bottom",
                         annotation_legend_side = annotationPosition,
                         column_title = title,
                         column_title_gp = grid::gpar(fontsize = 16))

}



#This function makes a HeatmapAnnotation object, and is a helper function for plotRegionsHeatmap and plotCNV

makeHeatmapAnnotation <- function(qseaSet,
                                  orientation,
                                  useGroupMeans = FALSE,
                                  sampleAnnotation = NULL){

  annotationColDf <- getAnnotation(qseaSet, sampleAnnotation = {{sampleAnnotation}}, useGroupMeans = useGroupMeans)

  if (is.null(annotationColDf)) {
    return(NULL)
  }

  annotationColDf <- annotationColDf %>%
    dplyr::mutate_if(is.character, as.factor)

  #Get all levels of all categorical variables and convert to color list
  levs <- annotationColDf %>%
    dplyr::select_if(is.factor) %>%
    purrr::map(function(x) levels(as.factor(x)))

  if (length(unlist(levs)) > 0) {


  col_list_cat <- levs %>%
    unlist() %>%
    length() %>%
    hues::iwanthue(plot=FALSE) %>%
    purrr::set_names(levs %>% unlist()) %>%
    utils::relist(levs) %>%
    purrr::map2(levs, purrr::set_names)
  } else {
    col_list_cat <- list()
  }

  #Make colour vectors for continuous variables
  annotationCol_numeric <- annotationColDf %>%
    dplyr::select_if(is.numeric)

  #Splitting into numeric variables with all non-negative values and those with negative values
  annotationCol_numeric_min_positive <- annotationCol_numeric %>%
    dplyr::select_if(function(x) min(x) >= 0)

  annotationCol_numeric_min_negative <- annotationCol_numeric %>%
    dplyr::select_if(function(x) min(x) < 0)

  colvecs_binary <- c("Reds","YlGnBu","YlOrBr","PuOr","Blues","Purples") %>%
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
    })

  annot <- ComplexHeatmap::HeatmapAnnotation(which = orientation,
                                             df    = annotationColDf,
                                             col   = annotationColors,
                                             annotation_legend_param = annotation_legend_param_ls,
                                             show_annotation_name    = FALSE)

  return(annot)
}

#' This function takes a qseaSet and a gene, and plots the expression across the gene as a heatmap
#' @param qseaSet The qseaSet object.
#' @param gene A gene to plot. Either a gene symbol or an ensembl ID.
#' @param normMethod Whether to plot nrpm values or beta values.
#' @param sampleAnnotation A data frame with annotations for the samples.
#' @param clusterNum A number of clusters to break the column dendrogram into.
#' @param clusterCols Whether to cluster the columns or not.
#' @param minDensity A minimum CpG density level to filter out windows with values lower than.
#' @param maxScale The maximum of the scale, not used when plotting beta values.
#' @param useGroupMeans Whether to use the group variable to collapse replicates.
#' @param upstreamDist Number of basepairs upstream of the gene to include.
#' @param downstreamDist Number of basepairs downstream of the gene to include.
#' @param minEnrichment Minimum enrichment factor for beta values, will give NAs below this.
#' @param scaleRows Whether to scale the rows of the heatmap
#' @param annotationColors A list with the colours to use for the column legend, to pass to pheatmap
#' @param mart A biomaRt mart object. If not supplied, will check the qseaSet, else will get a default for GRCh38/hg38 or hg19.
#' @param ... Additional arguments to pass to pheatmap.
#' @return A heatmap showing the methylation patterns across the gene of interest.
#' @export

plotGeneHeatmap <- function(qseaSet, gene, normMethod = "beta",
                            useGroupMeans = FALSE,
                            sampleAnnotation = NULL, minDensity = 0,
                            minEnrichment = 3, maxScale = 1, clusterNum = 2, annotationColors = NULL,
                            upstreamDist = 3000, scaleRows = FALSE, clusterCols = TRUE, mart = NULL,
                            downstreamDist = 1000, ...){

  annotationColDf = getAnnotation(qseaSet, sampleAnnotation = {{sampleAnnotation}}, useGroupMeans = useGroupMeans)

  if (!is.null(getMart(qseaSet))) { mart <- getMart(qseaSet) }

  if(is.null(mart) & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"Hsapiens") & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"hg38|GRCh38")){
    mart <- biomaRt::useMart('ensembl', dataset='hsapiens_gene_ensembl', host = "https://jul2022.archive.ensembl.org")
  } else if(is.null(mart) & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"Hsapiens") & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"hg19|GRCh37")) {
    mart <- biomaRt::useMart('ensembl', dataset='hsapiens_gene_ensembl', host = "https://feb2014.archive.ensembl.org")
  } else if(is.null(mart)){
    stop("Please specify a mart object for biomaRt.")
  }

  if (stringr::str_detect(gene,"^ENSG0")) {
    idType <- "ensembl_gene_id"
  } else{
    idType <- "hgnc_symbol"
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

  gene_details_gr <- gene_details %>%
    plyranges::as_granges()

  if (nrow(gene_details) != 1) {
    stop(glue::glue("Error: {nrow(gene_details)} genes matching this name found in {mart}."))
  }

  geneGR <- gene_details %>%
    plyranges::as_granges() %>%
    plyranges::anchor_start() %>%
    plyranges::stretch(upstreamDist) %>%
    plyranges::anchor_end() %>%
    plyranges::stretch(downstreamDist)

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

  if(useGroupMeans){
    colsToFind <- qseaSet %>% qsea::getSampleGroups() %>% names()
  } else {
    colsToFind <- qseaSet %>% qsea::getSampleNames()
  }

  numData <- dataTable %>%
    tibble::column_to_rownames("window") %>%
    dplyr::select(tidyselect::all_of(colsToFind)) %>%
    janitor::remove_empty(which = "cols", quiet = FALSE)

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
  colour_palette <- RColorBrewer::brewer.pal(name = "YlOrRd", n = 9)

  numData %>%
    as.matrix() %>%
    ComplexHeatmap::Heatmap(rect_gp = rectGpParam,
                            name = "Beta value",
                            cluster_rows = FALSE,
                            cluster_columns = clusterCols,
                            left_annotation = rowAnnot,
                            col = colour_palette,
                            show_row_names = F,
                            clustering_method_rows = "ward.D2",
                            clustering_method_columns = "ward.D2",
                            column_split = clusterNum,
                            column_title = NULL,
                            na_col = "lightgrey",
                            heatmap_legend_param = list(legend_direction = "vertical",
                                                        at = seq(0, maxScale, length.out = 6) %>% round(1))
                            ) %>%
    ComplexHeatmap::draw(heatmap_legend_side = "right",
                         annotation_legend_side = "right",
                         column_title = glue::glue("{stringr::str_to_title(normMethod)} values for {gene}"))

  return(invisible(numData))
}

#' This function makes a row annotation object for the gene heatmap, which hows information about gene regions and CpG density. It is a helper function for plotGeneHeatmap
#' @param rowAnnotationDF A dataframe of variables that provide information about the gene regions, which will be used to annotate windows within the heatmap

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

  colvecs_binary <- c("Reds","YlGnBu","YlOrBr","PuOr","Blues","Purples") %>%
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

#' This function takes a qseaSet and extracts the distribution of the genomic features of all the windows over a cutoff, and plots it
#' @param qseaSet The qseaSet object.
#' @param cutoff The cutoff to use on the windows for each sample
#' @param barType What type of bars to use for the plot (stack, dodge, fill)
#' @param normMethod Normalisation method to use
#' @return A plot of the distribution
#' @export
#'
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



#' This function takes a qseaSet and plots a correlation matrix
#' @param qseaSet The qseaSet object.
#' @param regionsToOverlap Regions to overlap
#' @param sampleAnnotation Columns of the sampleTable to use to annotate the heatmap
#' @param useGroupMeans Whether to use the replicate information or not (group column)
#' @param normMethod What normalisation method to use
#' @param annotationColors A list of custom colours to use for the annotation
#' @param minDensity Minimum CpG density to keep
#' @param minEnrichment Minimum number of reads for beta values
#' @param ... Other arguments for pheatmap
#' @return A table
#' @export
#'
plotCorrelationMatrix <- function(qseaSet, regionsToOverlap = NULL, useGroupMeans = FALSE, sampleAnnotation = NULL, normMethod = "nrpm",
                                  minEnrichment = 3, annotationColors = NA, minDensity = 0, ...){

  if (!is.null(regionsToOverlap)) {
    regionsToOverlap <- asValidGranges(regionsToOverlap)

    if (length(regionsToOverlap) == 0) {stop("No genomic regions given!")}

    qseaSet <- qseaSet %>%
      filterByOverlaps(regionsToOverlap = regionsToOverlap)
  }

  annotationDf = getAnnotation(qseaSet, sampleAnnotation = {{sampleAnnotation}}, useGroupMeans = useGroupMeans)

  dataTab <- qseaSet %>%
    filterWindows(CpG_density >= minDensity) %>%
    getDataTable(normMethod = normMethod,
                 useGroupMeans = useGroupMeans,
                 minEnrichment = minEnrichment
    ) %>%
    dplyr::mutate(window = paste0(seqnames, ":",start, "-",end)) %>%
    dplyr::mutate_all( ~ dplyr::case_when(!is.nan(.x) ~ .x)) # do something with NaN values?

  if(useGroupMeans){
    colsToFind <- qseaSet %>% qsea::getSampleGroups() %>% names()
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

#' This function takes a DMRtable as output by calculateDMRs (possibly after filtering), and generates an upset plot
#' @param DMRtable The data frame containing the differentially methylated regions.
#' @param string A string to subset the columns based on (always requires adjPval)
#' @param removeVS Whether to remove the string "_vs_" and everything after from the column names. For instance if all comparisons are against the same sample(s)
#' @param minAdjPval A minimum adjusted P value to consider the window as present
#' @param ... Other arguments to be passed to upset
#' @return An UpSet plot
#' @export
#'
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



#' Get an annotation data frame, for use in heatmap plotting functions
#'
#' Note that sampleAnnotation must be enclosed within double curly brackets when used.
#'
#' @param qseaSet A qseaSet
#' @param useGroupMeans Whether to use the "group" variable to collapse replicates
#' @param sampleAnnotation Columns of the sampleTable to use
#'
#' @return A data frame containing the annotation columns, ready for use in
getAnnotation <- function(qseaSet, useGroupMeans = FALSE, sampleAnnotation = NULL){

  if (rlang::quo_is_null(rlang::enquo(sampleAnnotation))) {
    return(NULL)
  }

  annotationColDf <- data.frame()

  if (!useGroupMeans) {
    annotationColDf <- qseaSet %>%
      qsea::getSampleTable() %>%
      dplyr::arrange(sample_name) %>%
      dplyr::select(!!!rlang::enquos(sampleAnnotation))

  } else if (useGroupMeans) {
    groupSampleTab <- qseaSet %>%
      qsea::getSampleTable() %>%
      dplyr::select(group, !!!rlang::enquos(sampleAnnotation))

    distinctGroupTab <- groupSampleTab %>%
      tibble::remove_rownames() %>%
      dplyr::distinct()

    if (nrow(distinctGroupTab) != (qseaSet %>% qsea::getSampleGroups() %>% length())) {
      problemGroups <- distinctGroupTab %>%
        dplyr::count(group) %>%
        filter(n > 1) %>%
        pull(group) %>%
        paste(., collapse = "; ")
      stop(glue::glue("Grouped annotation contains differing annotation in {problemGroups}"))
    }

    annotationColDf <- distinctGroupTab %>%
      arrange(group) %>%
      tibble::column_to_rownames("group") %>%
      as.data.frame()
  }

  if (ncol(annotationColDf) == 0) {
    annotationColDf <- NULL
  }

  return(annotationColDf)
}
