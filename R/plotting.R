
#' This function takes a qseaSet and a GRanges object, and plots the expression across the gene as a heatmap
#' @param qseaSet The qseaSet object.
#' @param signatureGR A genomic ranges to plot.
#' @param normMethod Whether to plot nrpm values or beta values.
#' @param annotationCol A data frame with annotations for the samples.
#' @param clusterNum A number of clusters to break the column dendrogram into.
#' @param maxScale The maximum of the scale, not used when plotting beta values.
#' @param clip A level to clip the data at, anything higher than this is replaced with the clip value.
#' @param minDensity A minimum CpG density level to filter out windows with values lower than.
#' @param clusterRows Whether to cluster the rows or not.
#' @param clusterCols Whether to cluster the columns or not.
#' @param clusterMethod What method to use to cluster the dendrograms.
#' @param description A string to include as part of the title.
#' @param annotationColors A list specifying some or all of the colours to use for the annotations.
#' @param minEnrichment Minimum enrichment factor for beta values, will give NAs below this.
#' @param ... Other arguments to pass to pheatmap.
#' @return A qseaSet object with the sampleTable enhanced with the information on number of reads etc
#' @export

plotGRangesHeatmap <- function(qseaSet, signatureGR, normMethod = "beta",
                               annotationCol = NULL,
                               annotationColors = NA, clusterRows = FALSE, clusterCols = TRUE,
                               minEnrichment = 3, maxScale = 5, clusterNum = 2, description = "",
                               clip = 1000000000, minDensity = 0,
                               clusterMethod = "ward.D2", ...){

  signatureGR <- plyranges::as_granges(signatureGR)

  if (length(signatureGR) == 0) {stop("No genomic regions given!")}

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

  remove_almost_empty_rows <- function(dat)  {
    mask_keep <- rowSums(is.na(dat)) != (ncol(dat) - 1)
    janitor:::remove_message(dat = dat, mask_keep = mask_keep, which = "rows", reason = "almost empty")
    return(dat[mask_keep, , drop = FALSE])
  }


  dataTab <- qseaSet %>%
    filterByOverlaps(windowsToKeep = signatureGR) %>%
    qsea::makeTable(., groupMeans = qsea::getSampleGroups(.), norm_methods = normMethod, minEnrichment = minEnrichment) %>%
    dplyr::filter(CpG_density >= minDensity) %>%
    dplyr::select(tidyselect::matches(normMethod)) %>%
    dplyr::select(-tidyselect::matches("PooledControl")) %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "_beta|_nrpm|_means")) %>%
    clipFn(a = 0, b = clip) %>%
    dplyr::mutate_all( ~ dplyr::case_when(!is.nan(.x) ~ .x)) %>%
    janitor::remove_empty(which = "cols", quiet = FALSE) %>%
    janitor::remove_empty(which = "rows", quiet = FALSE)

  if(ncol(dataTab) == 1) {
    clusterCols <- FALSE
  }

  if(clusterRows) {
    dataTab <- remove_almost_empty_rows(dataTab)
  }

  dataTab %>%
    pheatmap::pheatmap(annotation_col = annotationCol,
                       annotation_colors = annotationColors,
                       breaks = seq(0, maxScale, length.out = 10),
                       color = colour_palette,
                       cluster_rows = clusterRows,
                       cluster_cols = clusterCols,
                       show_rownames = FALSE,
                       treeheight_row = 0,
                       clustering_method = clusterMethod,
                       cutree_cols = clusterNum,
                       main = glue::glue("{normMethod} values for {description}."),
                       ...)


}

#' This function takes a qseaSet and a gene, and plots the expression across the gene as a heatmap
#' @param qseaSet The qseaSet object.
#' @param gene A gene to plot. Either a gene symbol or an ensembl ID.
#' @param normMethod Whether to plot nrpm values or beta values.
#' @param annotationCol A data frame with annotations for the samples.
#' @param clusterNum A number of clusters to break the column dendrogram into.
#' @param clusterCols Whether to cluster the columns or not.
#' @param minDensity A minimum CpG density level to filter out windows with values lower than.
#' @param maxScale The maximum of the scale, not used when plotting beta values.
#' @param upstreamDist Number of basepairs upstream of the gene to include.
#' @param downstreamDist Number of basepairs downstream of the gene to include.
#' @param minEnrichment Minimum enrichment factor for beta values, will give NAs below this.
#' @param scaleRows Whether to scale the rows of the heatmap
#' @param annotationColors A list with the colours to use for the column legend, to pass to pheatmap
#' @param mart A biomaRt mart object. If not supplied, will check the qseaSet, else will get a default for GRCh38/hg38 or hg19.
#' @param ... Additional arguments to pass to pheatmap.
#' @return A qseaSet object with the sampleTable enhanced with the information on number of reads etc
#' @export

plotGeneHeatmap <- function(qseaSet, gene, normMethod = "beta",
                             annotationCol = NULL, minDensity = 0,
                             minEnrichment = 3, maxScale = 1, clusterNum = 2, annotationColors = NULL,
                             upstreamDist = 3000, scaleRows = FALSE, clusterCols = TRUE, mart = NULL,
                             downstreamDist = 1000, ...){

  if (!is.null(annotationCol) & !is.data.frame(annotationCol) & all(is.character(annotationCol))) {
    annotationCol <- getAnnotationDataFrame(qseaSet, !!!rlang::syms(annotationCol))
  } else if (is.null(annotationCol) ){
    annotationCol <- NA
  }

  if(!is.null(getMart(qseaSet))){ mart <- getMart(qseaSet) }

  if(is.null(mart) & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"Hsapiens") & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"hg38|GRCh38")){
    mart <- biomaRt::useMart('ensembl', dataset='hsapiens_gene_ensembl')
  } else if(is.null(mart) & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"Hsapiens") & stringr::str_detect(qseaSet %>% qsea:::getGenome(),"hg19|GRCh37")) {
    mart <- biomaRt::useMart('ensembl', dataset='hsapiens_gene_ensembl', host = "https://feb2014.archive.ensembl.org")
  } else if(is.null(mart)){
    stop("Please specify a mart object for biomaRt.")
  }

  if(stringr::str_detect(gene,"^ENSG0")){
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

  qseaSetChr <- qseaSet %>% qsea::getRegions() %>% GenomeInfoDb::seqinfo() %>% GenomeInfoDb::seqnames() %>% stringr::str_detect("chr") %>% any()
  windowsChr <- gene_details %>% pull(seqnames) %>% stringr::str_detect("chr") %>% any()

  if(qseaSetChr & !windowsChr){
    gene_details <- gene_details %>%
      dplyr::mutate(seqnames = paste0("chr", seqnames))
  }

  gene_details_gr <- gene_details %>%
    plyranges::as_granges()

  if(nrow(gene_details) != 1){
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
    filterByOverlaps(windowsToKeep = geneGR) %>%
    qsea::makeTable(., groupMeans = qsea::getSampleGroups(.),
                    norm_methods = normMethod, minEnrichment = minEnrichment) %>%
    dplyr::filter(CpG_density >= minDensity) %>%
    dplyr::mutate(window = paste0(chr, ":", window_start, "-", window_end))

  geneStrand <- gene_details_gr %>% tibble::as_tibble() %>% dplyr::pull(strand)

  annoRow <- dataTable %>%
    dplyr::rename(seqnames = chr, start = window_start, end = window_end) %>%
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

  if (scaleRows) {
    scaleRows <- "row"
    colour_palette <- RColorBrewer::brewer.pal(name = "RdBu", n = 9)
    breaks <- seq(-3, 3, length.out = 10)
  } else{
    scaleRows <- "none"
    colour_palette <- RColorBrewer::brewer.pal(name = "YlOrRd", n = 9)
    breaks <- seq(0, maxScale, length.out = 10)
  }

  dat <- dataTable %>%
    dplyr::select(window, dplyr::matches(normMethod)) %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "_beta|_nrpm|_means")) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("window") %>%
    janitor::remove_empty(which = "cols", quiet = FALSE)

  # remove_almost_empty_cols <- function(dat)  {
  #   mask_keep <- colSums(is.na(dat)) != (nrow(dat) - 1)
  #   janitor:::remove_message(dat = dat, mask_keep = mask_keep, which = "cols", reason = "almost empty")
  #   return(dat[, mask_keep, drop = FALSE])
  # }
  #
  # if(clusterCols) {
  #   dat <- remove_almost_empty_cols(dat)
  # }
  #

  dat %>%
    pheatmap::pheatmap(annotation_col = annotationCol,
                       annotation_row = annoRow,
                       annotation_colors = annotationColors,
                       breaks = breaks,
                       color = colour_palette,
                       cluster_rows = FALSE,
                       cluster_cols = clusterCols,
                       scale = scaleRows,
                       show_rownames = FALSE,
                       treeheight_row = 0,
                       clustering_method = "ward.D2",
                       cutree_cols = clusterNum,
                       main = glue::glue("{str_to_title(normMethod)} values for {gene}."),
                       ...
    )

  return(invisible(dat))

}




#' This function takes a qseaSet and plots a PCA
#' @param qseaSet The qseaSet object.
#' @param signatureGR A GRanges (or data frame coercible to one) with a subset of windows to calculate the PCA on.
#' @param plotComponents Which of the PCA components to plot default c(1,2)
#' @param plotColour Column of the sampleTable to use to set the colour of the points
#' @param plotShape Column of the sampleTable to use to set the shape of the points
#' @param batchVariable A variable to use for batch correction using limma. Can be overaggressive and remove all variation. Be careful with this!
#' @param showNames Whether to show the names of the points
#' @param normMethod Normalisation method to use, nrpm or beta.
#' @param minRowSum A minimum number of reads to be included in the PCA
#' @param minEnrichment Minimal number of expected reads for a fully methylated window (for transformation to absolute methylation level)
#' @param topVar Take only the n most variable windows in the PCA. If NULL, all valid windows.
#' @param plotTitle Title for the plot
#' @param returnDataOnly Return just the data frame behind the PCA
#' @return A plot of the PCA
#' @export
#'
plotQseaPCA <- function(qseaSet,
                        signatureGR = NULL,
                        plotComponents = c(1,2),
                        plotColour = "type",
                        plotShape = "experiment",
                        batchVariable = NULL,
                        showNames = FALSE,
                        normMethod = "beta",
                        minRowSum = 20,
                        minEnrichment = 3,
                        topVar = NULL,
                        plotTitle = "",
                        returnDataOnly = FALSE
){

  qseaSet <- dropPooledControl(qseaSet)

  if (!is.null(signatureGR)) {
    qseaSet <- filterByOverlaps(qseaSet, signatureGR)
  }

  pc1 <- plotComponents[1]
  pc2 <- plotComponents[2]


  pca <- getPCAwithBatch(qseaSet,
                         normMethod = normMethod,
                         topVar = topVar,
                         batchVariable = batchVariable,
                         minRowSum = minRowSum,
                         minEnrichment = minEnrichment)

  svd <- pca@svd
  propVar <- (svd$d^2)/sum(svd$d)
  propVar <- round(propVar/sum(propVar)*100,2)

  plotData <- tibble::tibble(
    x = svd$v[,pc1]*svd$d[pc1],
    y = svd$v[,pc2]*svd$d[pc2],
    sample_name = qsea::getSampleNames(qseaSet)
  ) %>%
    dplyr::left_join(qsea::getSampleTable(qseaSet), by = "sample_name")

  batchTitleString <- ifelse(!is.null(batchVariable),"Batch-corrected ","")

  if(returnDataOnly){return(plotData)}

  p <- ggplot2::ggplot(plotData, ggplot2::aes_string("x", "y",
                                                label = "sample_name",
                                                colour = plotColour,
                                                shape = plotShape)) +
    ggplot2::geom_point() +
    ggplot2::xlab(glue::glue("PC{pc1} ({propVar[pc1]}%)")) +
    ggplot2::ylab(glue::glue("PC{pc2} ({propVar[pc2]}%)")) #+
    #ggplot2::ggtitle(glue::glue("{plotTitle}"),
    #                 subtitle = glue::glue("{batchTitleString}PCA on {length(pca@factor_names)} windows, using {normMethod} values."))

  if (showNames) {

      if (!requireNamespace("ggrepel", quietly = TRUE)) {
        message("Package \"ggrepel\" is recommended to repel labels. Using default method.")
        p <- p + ggplot2::geom_text()
      }
      p <- p + ggrepel::geom_text_repel()
    }

  p
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
#' @param normMethod What normalisation method to use
#' @param annotationColors A data frame for use of annotation of the rows and columns in pheatmap
#' @param ... Other arguments for pheatmap
#' @return A table
#' @export
#'
plotCorrelationMatrix <- function(qseaSet, ..., normMethod = "nrpm", annotationColors = NA){

  annotationCol = getAnnotationDataFrameIndividual(qseaSet, ...)

  qseaSet %>%
    dropPooledControl() %>%
    qsea::makeTable(samples = qsea::getSampleNames(.), norm_methods = normMethod) %>%
    dplyr::select(tidyselect::matches(normMethod)) %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "_beta|_nrpm")) %>%
    stats::cor(use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    pheatmap::pheatmap(display_numbers = TRUE,
                       color = RColorBrewer::brewer.pal(name = "YlOrRd", n = 9),
                       annotation_row = annotationCol,
                       annotation_col = annotationCol,
                       annotation_colors = annotationColors)
}

#' This function takes a DMRtable as output by calculateDMRs (possibly after any filtering), and generates an upset plot
#' @param DMRtable The data frame containing the differentially methylated regions.
#' @param string A string to subset the columns based on (always requires adjPval)
#' @param removeVS Whether to remove the string "_vs_" and everything after from the column names. For instance if all comparisons are against the same sample(s)
#' @param minAdjPval A minimum adjusted P value to consider the window as present
#' @param ... Other arguments to be passed to upset
#' @return A table
#' @export
#'
plotDMRUpSet <- function(DMRtable, string = NULL, removeVS = TRUE, minAdjPval = 0.05, ...){

  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop(
      "Package \"UpSetR\" must be installed to use this function.",
      call. = FALSE
    )
  }

  temp <- DMRtable %>%
    dplyr::select(tidyselect::matches("adjPval")) %>%
    {if (removeVS) {dplyr::rename_with(., ~stringr::str_remove(.,"_vs_.*")) } else .} %>%
    {if (!is.null(string)) {dplyr::select(., tidyselect::matches(string))} else . }

  if (ncol(temp) < 2) {stop(glue::glue("Only {ncol(temp)} column remaining"))}

  purrr::map(stats::setNames(colnames(temp),colnames(temp)),function(x){  which(temp[,x,drop = FALSE] <= minAdjPval)}) %>%
    UpSetR::fromList() %>%
    UpSetR::upset(nsets = ncol(.), nintersects = 100, order.by = "freq", text.scale	 = 1.8, ...)
}
