
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
    dplyr::mutate_all( ~ case_when(!is.nan(.x) ~ .x)) %>%
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
#' @param gene A gene to plot. Currently as a gene id only, in a dubious way.
#' @param normMethod Whether to plot nrpm values or beta values.
#' @param annotationCol A data frame with annotations for the samples.
#' @param clusterNum A number of clusters to break the column dendrogram into.
#' @param clusterCols Whether to cluster the columns or not.
#' @param maxScale The maximum of the scale, not used when plotting beta values.
#' @param upstreamDist Number of basepairs upstream of the gene to include.
#' @param downstreamDist Number of basepairs downstream of the gene to include.
#' @param minEnrichment Minimum enrichment factor for beta values, will give NAs below this.
#' @param scaleRows Whether to scale the rows of the heatmap
#' @param annotationColors A list with the colours to use for the column legend, to pass to pheatmap
#' @return A qseaSet object with the sampleTable enhanced with the information on number of reads etc
#' @export

plotGeneHeatmap <- function(qseaSet, gene, normMethod = "beta", annotationCol = NULL,
                            minEnrichment = 3, maxScale = 1, clusterNum = 2, annotationColors = NA,
                            upstreamDist = 3000, scaleRows = FALSE, clusterCols = TRUE,
                            downstreamDist = 1000){

  if (is.null(annotationCol) & ("type" %in% (qseaSet %>% qsea::getSampleTable() %>% colnames()))) {
    annotationCol <- getAnnotationDataFrame(qseaSet, type)
  }

  # if (typeof(annotationCol) == "character") {
  #   annotationCol <- getAnnotationDataFrame(qseaSet, annotationCol)
  # }


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



  geneGR <- mesa::GRCh38.99.gtf %>%
    dplyr::filter(type == "gene", gene_name == gene) %>%
    plyranges::anchor_start() %>%
    plyranges::stretch(upstreamDist) %>%
    plyranges::anchor_end() %>%
    plyranges::stretch(downstreamDist)

  if (length(geneGR) == 0) {stop("No genomic region found!")}
  if (length(geneGR) > 1) {stop("Multiple genomic regions found!")}

  if (normMethod == "beta") {maxScale = 1}

  # fivePrimes <- GRCh38.99.gtf %>%
  #   filter(gene_name == gene, type == "five_prime_utr",  tag == "basic")
  #
  # threePrimes <- GRCh38.99.gtf %>%
  #   filter(gene_name == gene, type == "three_prime_utr",  tag == "basic")
  #
  # exons <- GRCh38.99.gtf %>%
  #   filter(gene_name == gene, type == "exon",  tag == "basic")
  #
  # genebody <- GRCh38.99.gtf %>%
  #   filter(gene_name == gene, type == "gene",  tag == "basic")
  #
  # annotatedData2 <- qseaSet %>%
  #   filterByOverlaps(windowsToKeep = geneGR) %>%
  #   qsea::makeTable(., groupMeans = qsea::getSampleGroups(.),
  #                   norm_methods = normMethod, minEnrichment = minEnrichment) %>%
  #   mutate(seqnames = chr, start = window_start, end = window_end) %>%
  #   plyranges::as_granges() %>%
  #   plyranges::mutate(fivePrime = plyranges::count_overlaps(.,fivePrimes),
  #                     threePrime = plyranges::count_overlaps(.,threePrimes),
  #                     exon = plyranges::count_overlaps(.,exons),
  #                     genebody = plyranges::count_overlaps(.,genebody)) %>%
  #   as_tibble() %>%
  #   dplyr::mutate(position = if_else(genebody > 0 , "intron", ""),
  #          position = if_else(exon > 0, "exon", position),
  #          position = if_else(threePrime > 0, "3'UTR", position),
  #          position = if_else(fivePrime > 0, "5'UTR", position)
  #         )

  annotatedData <- qseaSet %>%
    filterByOverlaps(windowsToKeep = geneGR) %>%
    qsea::makeTable(., groupMeans = qsea::getSampleGroups(.),
                    norm_methods = normMethod, minEnrichment = minEnrichment) %>%
    annotateData() %>%
    dplyr::mutate(window = paste0(chr, ":", window_start, "-", window_end))

  annoRow <- annotatedData %>%
    #dplyr::select(window, CpG_density) %>%
    dplyr::select(window, CpG_density, annotation) %>%
    dplyr::mutate(annotation = stringr::str_replace(annotation, "Exon .*", "Exon")) %>%
    dplyr::mutate(annotation = stringr::str_replace(annotation, "Intron .*", "Intron")) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("window") %>%
    as.data.frame()

  if (scaleRows) {
    scaleRows <- "row"
    colour_palette <- RColorBrewer::brewer.pal(name = "RdBu", n = 9)
    breaks <- seq(-3, 3, length.out = 10)
  } else{
    scaleRows <- "none"
    colour_palette <- RColorBrewer::brewer.pal(name = "YlOrRd", n = 9)
    breaks <- seq(0, maxScale, length.out = 10)
  }

  annotatedData %>%
    dplyr::select(window, dplyr::matches(normMethod)) %>%
    dplyr::select(-dplyr::matches("PooledControl")) %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "_beta|_nrpm|_means")) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("window") %>%
    janitor::remove_empty(which = "cols", quiet = FALSE) %>%
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
                       main = glue::glue("Heatmap of {normMethod} values for {gene}."))

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
