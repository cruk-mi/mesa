#' Calculate and add CNV based on HMMcopy, on individual samples
#'
#' This calculates and adds CNV to a qseaSet based on just using HMMcopy on each sample, with default parameters.
#' Only works for hg38/GRCh38!
#'
#' @param qs The qseaSet object
#' @param inputColumn Which column of the sample table contains the link to the bam file
#' @param windowSize What window size to use (only 50000, 500000 or 1000000 work at the moment)
#' @param parallel Whether to use parallelisation
#' @param fragmentLength Fragment length to extend reads to if not paired
#' @param minMapQual Minimum mapping quality to include a read (on either end if proper pair)
#' @param minInsertSize For paired reads, only keep them if they are above a minimum length. Can be used for cfDNA size selection. Applies to Input samples as well as MeCap.
#' @param maxInsertSize For paired reads, only keep them if they are below a maximum length. Can be used for cfDNA size selection. Applies to Input samples as well as MeCap.
#' @param properPairsOnly Whether to only keep properly paired reads, or to keep high-quality (MAPQ 30+) unpaired R1s as well. Set to TRUE for size selection.
#' @param minReferenceLength A minimum distance on the genome to keep the read. bwa by default gives 19bp as minimum for a read, which is quite short.
#' @param plotDir A directory to export individual HMMcopy plots to.
#' @return A qseaGLM object
#' @export

addHMMcopyCNV <- function(qs, inputColumn = "input_file", windowSize = 1000000, fragmentLength = NULL,
                          plotDir = NULL, parallel = TRUE,
                          maxInsertSize = 1000,
                          minInsertSize = 50,
                          minReferenceLength = 30,
                          minMapQual = 30,
                          properPairsOnly = FALSE
                          )
{

  #TODO something with zygosity

  if(windowSize == 1000000) {

    gcGR <- gc_hg38_1000kb
    mapGR <- map_hg38_1000kb

  } else if (windowSize == 500000) {

    gcGR <- gc_hg38_500kb
    mapGR <- map_hg38_500kb

  } else if (windowSize == 50000) {

    gcGR <- gc_hg38_50kb
    mapGR <- map_hg38_50kb

  } else if (windowSize == 10000) {

    gcGR <- gc_hg38_10kb
    mapGR <- map_hg38_10kb

  } else {stop(glue::glue("Window size must be one of 10000, 50000, 500000 or 1000000 currently."))}

  CNV_Regions = qsea:::makeGenomeWindows(qsea:::getGenome(qs), as.character(qsea::getChrNames(qs)), windowSize)
  CpGpos = GenomicRanges::GRanges(seqinfo = GenomicRanges::seqinfo(CNV_Regions))

  if (parallel) {
    BPPARAM = BiocParallel::bpparam()
    message("Scanning up to ", BiocParallel::bpnworkers(BPPARAM), " files in parallel")
  }  else {BPPARAM = BiocParallel::SerialParam()}

  bamOutList <- BiocParallel::bplapply(X = qsea::getSampleTable(qs) %>% dplyr::pull(inputColumn),
                                       FUN = getBamCoveragePairedAndUnpairedR1, regions = CNV_Regions,
                                       fragmentLength = fragmentLength, maxInsertSize = maxInsertSize,
                                       minInsertSize = minInsertSize, minReferenceLength = minReferenceLength,
                                       minMapQual = minMapQual, properPairsOnly = properPairsOnly, BPPARAM = BPPARAM)

  names(bamOutList) <- qsea::getSampleTable(qs)$sample_name

  libraries = bamOutList %>%
    purrr::map_dfr(~ purrr::pluck(.,"library")) %>%
    as.data.frame()

  rownames(libraries) <- qsea::getSampleTable(qs)$sample_name

  qs = qsea:::setLibrary(qs, "input_file", libraries)

  counts =  bamOutList %>%
    purrr::map_dfc(function(x) {
      purrr::pluck(x,"regionsGR") %>%
        tibble::as_tibble() %>%
        dplyr::pull(counts)  }
    ) %>% as.matrix()
  colnames(counts) <- qsea::getSampleTable(qs)$sample_name

  nf = apply(X = counts, MARGIN = 2, FUN = stats::median, na.rm = TRUE)
  if (any(nf < 100, na.rm = TRUE)) {
    warning("Low coverage in CNV files. Consider larger CNV_window_size")
  }

  CNV_RegionsWithReads <- CNV_Regions %>%
    tibble::as_tibble() %>%
    dplyr::rename(chr = seqnames) %>%
    dplyr::bind_cols(reads = tibble::as_tibble(counts)) %>%
    dplyr::left_join(gcGR, by = c("chr", "start", "end")) %>%
    dplyr::left_join(mapGR, by = c("chr", "start", "end")) %>%
    dplyr::mutate(chr = as.factor(chr)) %>%
    data.table::data.table()

  #try saving the raw CNV data in the qseaSet.
  #TODO would be better as a slot of its own, but can't seem to manage that into the S4 class
  qs@libraries$input_counts <- CNV_RegionsWithReads

  mergedCNV <- purrr::map(qsea::getSampleTable(qs)$sample_name, ~ runHMMCopy(CNV_RegionsWithReads, ., plotDir = plotDir)) %>%
    purrr::reduce(plyranges::join_overlap_inner)

  qs = qsea:::setCNV(qs, mergedCNV)

  if (length(qsea::getOffset(qs)) > 0 && !all(is.na(qsea::getOffset(qs))))
    warning("Consider recalculating offset based on new CNV values")
  if ("seqPref" %in% names(GenomicRanges::mcols(qsea::getRegions(qs))))
    warning("Consider recalculating sequence ", "preference based on new CNV values")
  return(qs)
}


#' Function to call HMMcopy on the individual reads for a sample
#'
#' This calculates and adds CNV to a qseaSet based on just using HMMcopy on each sample, with default parameters.
#' Only works for hg38/GRCh38!
#'
#' @param CNV_RegionsWithReads Granges object with reads columns
#' @param colname Which column of the sample to use
#' @param plotDir A directory to export individual HMMcopy plots to.
#' @return A GRanges object with the result of calling HMMcopy.
#' @export
runHMMCopy <- function(CNV_RegionsWithReads, colname, plotDir = NULL){

  if(!is.null(plotDir)){
      dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)
  }

  correctOutput <- CNV_RegionsWithReads %>%
    dplyr::select(chr, start, end, width, strand, gc, map, tidyselect::any_of(colname)) %>%
    dplyr::rename(reads = colname) %>%
    HMMcopy::correctReadcount(mappability = 0.9, samplesize = 50000, verbose = FALSE)

  initParam <- HMMcopy::HMMsegment(correctOutput, param = NULL, autosomes = NULL,
                          maxiter = 50, getparam = TRUE, verbose = FALSE)

  hmmseg <- HMMcopy::HMMsegment(correctOutput, param = initParam, autosomes = NULL,
                       maxiter = 50, verbose = FALSE)

  endMusDf <- round(hmmseg$mus[,ncol(hmmseg$mus)], 2) %>%
    tibble::enframe(name = "state", value = "CNV") %>%
    dplyr::mutate(state = as.character(state))

  out <- hmmseg$segs %>%
    dplyr::left_join(endMusDf, by = "state") %>%
    dplyr::mutate(chr = factor(chr, levels = 1:22)) %>%
    dplyr::arrange(chr, start) %>%
    dplyr::rename(seqnames = chr) %>%
    dplyr::select(seqnames, start, end, CNV) %>%
    plyranges::as_granges()

  colnames(GenomicRanges::mcols(out)) <- stringr::str_replace(colnames(GenomicRanges::mcols(out)),"CNV",colname)

  tttt <- hmmseg$state %>%
    tibble::as_tibble() %>%
    dplyr::mutate(state = as.character(value)) %>%
    dplyr::left_join(endMusDf, by = "state")

  if (!is.null(plotDir)) {

    grDevices::pdf(glue::glue("{plotDir}/{colname}.pdf"), width = 9, height = 9)
    plot(correctOutput$copy, col = "black",
         xlab = "Window (over chr1 to chr22)",
         ylab = "Log 2 Ratio compared to 2 chromosomes")
    graphics::points(tttt$CNV,col = "red", lwd = 0.5)
    #points(as.data.frame(qseaM023@cnv)[,colname], col = "blue", lwd = 0.5)
    #points(as.data.frame(qseaM023recalInput@cnv)[,colname], col = "green", lwd = 0.5)
    graphics::title(main = glue::glue("Copy number estimated via HMMcopy."),
          sub = "Black circles are HMMcopy mappability and gc normalised values.")
    grDevices::dev.off()

  }

  CNV_RegionsWithReads %>%
    dplyr::select(chr, start, end) %>%
    dplyr::rename(seqnames = chr) %>%
    plyranges::as_granges() %>%
    plyranges::join_overlap_inner(out) %>%
    return()

}

#' This function takes a qseaSet and plots a heatmap of the calculated CNV
#' @param qseaSet The qseaSet object.
#' @param annotationCol A data frame with annotations for the samples (can be made with getAnnotationDataFrameIndividual)
#' @param annotationColors A list with the colours to use for the column legend, to pass to pheatmap
#' @param clusterRows Whether to cluster the rows of the heatmap
#' @return A heatmap with the calculated number of chromosomes for each samples
#' @export
#'
plotCNVheatmap <- function(qseaSet,
                           ...,
                           annotationCol = NULL,
                           annotationColors = NA,
                           clusterRows = TRUE){

  if(is.null(annotationCol)){
    annotationCol = getAnnotationDataFrameIndividual(qseaSet, ...)
  }

  # if the annotation is empty then set to NA else pheatmap complains
  if(ncol(annotationCol) == 0){
    annotationCol <- NA
  }

  if (is.list(annotationColors)) {

    namesVec <- names(annotationColors)
    names(namesVec) <- namesVec

    annotationColors <- purrr::map(namesVec,
      function(x){
        usedValues <- qseaSet %>%
          qsea::getSampleTable() %>%
          dplyr::pull(x) %>%
          unique()
        return(x = annotationColors[[x]][usedValues])
      }
    )

    if(!("chr" %in% names(annotationColors))){
      annotationColors$chr <- c("1" = "black", "0" = "darkgrey")
    }

  }



  # Take the CNV and generate a data frame with the annotation
  rowAnno <- qseaSet %>%
    qsea::getCNV() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = factor(seqnames,
                                    levels = sort(as.numeric(levels(seqnames))))) %>%
    dplyr::mutate(window = paste0(seqnames, ":",start, "-",end)) %>%
    #dplyr::mutate(seqnames = factor(seqnames, levels = c(1:22))) %>%
    dplyr::mutate(chr = as.numeric(seqnames) %% 2 ) %>%
    dplyr::mutate(chr = factor(chr,levels = c(1,0))) %>%
    dplyr::arrange(seqnames) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("window") %>%
    dplyr::select(chr)

  CNVmatrix <- qseaSet %>%
    qsea::getCNV() %>%
    tibble::as_tibble() %>%
    dplyr::select(-tidyselect::matches("PooledControl")) %>%
    dplyr::mutate(window = paste0(seqnames, ":",start, "-",end)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("window") %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>%
    as.matrix() %>%
    t()
    #{2*(2^.)} %>%

  minmax <- max(max(CNVmatrix),abs(min(CNVmatrix)),1)

  CNVmatrix %>%
    ComplexHeatmap::pheatmap(cluster_rows = clusterRows,
                       cluster_cols = FALSE,
                       show_colnames = FALSE,
                       annotation_row = annotationCol,
                       annotation_colors = annotationColors,
                       #breaks = 2*(2^seq(from = 0, to = 4, length.out = 12)),
                       breaks = seq(from = -minmax, to = minmax, length.out = 12),
                       #labels_row = stringr::str_remove(rownames(.), "_.*"),
                       annotation_col = rowAnno,
                       gaps_col = rowAnno %>% dplyr::pull(chr) %>% diff() %>% {which(. != 0)},
                       color = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11))
    )

}

#' This function takes a qseaSet and removes the CNV data.
#' @param qseaSet The qseaSet object.
#' @return A qseaSet with the CNV data removed
#' @export
#'
removeCNV <- function(qseaSet){

  qseaSet@cnv <- qseaSet@cnv %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(-c(seqnames, start, end, width, strand),~.*0)) %>%
    plyranges::as_granges()

  return(qseaSet)

}
