#' Add per-sample CNV calls to a qseaSet using HMMcopy
#'
#' Runs a CNV pipeline on each sample: counts fragments in fixed-size genomic
#' windows, normalises for GC and mappability (HMMcopy), segments, and stores
#' per-window copy-number in the \code{qseaSet}. Supports parallel execution.
#' Only GRCh38/hg38 is currently supported.
#'
#' @details
#' For each sample, coverage is computed with
#' \code{\link{getBamCoveragePairedAndUnpairedR1}} over genome windows generated
#' by \code{qsea:::makeGenomeWindows()}, normalised with HMMcopy
#' (\code{HMMcopy::correctReadcount()} + \code{HMMcopy::HMMsegment()}), and the
#' resulting CNV tracks are merged across samples into the qseaSet.
#'
#' @param qs A \code{qseaSet} object.
#' @param inputColumn Character(1). Column name in the sample table with BAM paths (e.g. `"input_file"`).
#' @param windowSize Integer(1). Window size in bp; allowed values: 50000, 500000, 1000000.
#' @param parallel Logical(1). Use \pkg{BiocParallel} if a parallel backend is registered.
#' @param fragmentLength Integer(1) or \code{NULL}. Extend unpaired R1 to this length; if \code{NULL},
#'   estimate from proper pairs.
#' @param minMapQual Integer(1). Minimum MAPQ for either end (pairs) or R1.
#' @param minInsertSize,maxInsertSize Integer(1). Absolute insert-size bounds for proper pairs.
#' @param properPairsOnly Logical(1). If \code{TRUE}, ignore unpaired R1; if \code{FALSE}, combine pairs + high‑quality R1.
#' @param minReferenceLength Integer(1). Minimum reference span for R1 (unpaired).
#' @param plotDir Character(1) or \code{NULL}. If given, PDF plots per sample are written here.
#' @param hmmCopyGC,hmmCopyMap \linkS4class{GRanges}. Precomputed GC and mappability tracks binned at \code{windowSize}.
#'
#' @return The input \code{qseaSet} with:
#' \itemize{
#'   \item CNV per window stored via \code{qsea:::setCNV()};
#'   \item \code{@libraries$input_counts} containing raw per-window counts + GC/map;
#'   \item library table updated via \code{qsea:::setLibrary()}.
#' }
#'
#' @seealso \code{\link{runHMMCopy}}, \code{\link{plotCNVheatmap}},
#'   \code{\link{getBamCoveragePairedAndUnpairedR1}}, \pkg{HMMcopy}
#' @family CNV
#' 
#' @examples
#' # Example workflow (non-running here; requires BAMs and GC/map tracks)
#' # data(exampleTumourNormal, package = "mesa")
#' # qs <- exampleTumourNormal
#' #
#' # # Provide GC and mappability tracks at the desired window size
#' # # (e.g. precomputed 1 Mb hg38 tracks as GRanges with columns 'gc' and 'map')
#' # gc_track  <- hmmcopy_gc_1Mb_hg38   # user-supplied GRanges
#' # map_track <- hmmcopy_map_1Mb_hg38  # user-supplied GRanges
#' #
#' # # Make sure the sample table has a BAM column (e.g., 'input_file')
#' # # st <- qsea::getSampleTable(qs); stopifnot("input_file" %in% names(st))
#' #
#' # BiocParallel::register(BiocParallel::SerialParam())
#' # qs2 <- addHMMcopyCNV(
#' #   qs, inputColumn = "input_file", windowSize = 1e6,
#' #   hmmCopyGC = gc_track, hmmCopyMap = map_track,
#' #   properPairsOnly = TRUE, parallel = FALSE, plotDir = tempdir()
#' # )
#' # # Per-window CNV now available:
#' # # qsea::getCNV(qs2)
#' 
#' @export
addHMMcopyCNV <- function(qs, inputColumn = "input_file", windowSize = 1000000, fragmentLength = NULL,
                          plotDir = NULL,
                          parallel = getMesaParallel(),
                          maxInsertSize = 1000,
                          minInsertSize = 50,
                          minReferenceLength = 30,
                          minMapQual = 30,
                          properPairsOnly = FALSE,
                          hmmCopyGC = NULL,
                          hmmCopyMap = NULL
                          )
{

  if (is.null( hmmCopyGC)) {
    stop("GC track for hmmcopy must be provided.")
  }

  if (is.null( hmmCopyMap)) {
    stop("Mapability track for hmmcopy must be provided.")
  }

  hmmCopyGC <- asValidGranges(hmmCopyGC)
  hmmCopyMap <- asValidGranges(hmmCopyMap)

  #TODO something with zygosity

  CNV_Regions = qsea:::makeGenomeWindows(qsea:::getGenome(qs), as.character(qsea::getChrNames(qs)), windowSize)
  CpGpos = GenomicRanges::GRanges(seqinfo = GenomicRanges::seqinfo(CNV_Regions))

  if (parallel) {
    BPPARAM = BiocParallel::bpparam()
    message("Scanning up to ", BiocParallel::bpnworkers(BPPARAM), " files in parallel")
  }  else {BPPARAM = BiocParallel::SerialParam()}

  bamOutList <- BiocParallel::bplapply(X = qsea::getSampleTable(qs) %>% dplyr::pull(inputColumn),
                                       FUN = getBamCoveragePairedAndUnpairedR1,
                                       BSgenome = qsea:::getGenome(qs),
                                       regions = CNV_Regions,
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

  GenomicRanges::mcols(CNV_Regions) <- tibble::as_tibble(counts)

  CNV_RegionsWithReads <- CNV_Regions %>%
    plyranges::join_overlap_left(hmmCopyGC) %>%
    plyranges::join_overlap_left(hmmCopyMap)  %>%
    tibble::as_tibble() %>%
    dplyr::rename(chr = seqnames) %>%
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


#' Run HMMcopy on per-window reads for a single sample
#'
#' Given a table of per-window counts with GC and mappability, run HMMcopy
#' correction and segmentation, and return a GRanges of CNV calls for the
#' selected sample/column.
#'
#' @param CNV_RegionsWithReads A \code{data.frame} / \code{data.table} / \linkS4class{GRanges}
#'   coercible object with columns:
#'   \code{chr}, \code{start}, \code{end}, \code{width}, \code{strand}, \code{gc}, \code{map},
#'   and one sample-specific \code{reads} column (passed by name via \code{colname}).
#' @param colname Character(1). The sample column to use as read counts.
#' @param plotDir Character(1) or \code{NULL}. If provided, writes a PDF plot named \code{<colname>.pdf}.
#'
#' @return A \linkS4class{GRanges} with columns named \code{<colname>} containing
#'   per-window CNV estimates (HMM state means on the log2 scale).
#'
#' @seealso \code{\link{addHMMcopyCNV}}, \pkg{HMMcopy}
#' @family CNV
#'
#' @examples
#' # Minimal skeleton (non-running): build the required table, then call:
#' # tbl <- data.frame(
#' #   chr = rep(paste0("chr", 1:2), each = 5),
#' #   start = seq(1, by = 1e6, length.out = 10),
#' #   end   = seq(1e6, by = 1e6, length.out = 10),
#' #   width = 1e6, strand = "*",
#' #   gc = runif(10, 0.3, 0.7), map = runif(10, 0.8, 1.0),
#' #   sampleA = rpois(10, 500)
#' # )
#' # gr_cnv <- runHMMCopy(tbl, colname = "sampleA", plotDir = tempdir())
#' # gr_cnv
#' 
#' @export
runHMMCopy <- function(CNV_RegionsWithReads, colname, plotDir = NULL){

  if(!is.null(plotDir)){
      dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)
  }

  correctOutput <- CNV_RegionsWithReads %>%
    dplyr::select(chr, start, end, width, strand, gc, map, tidyselect::any_of(colname)) %>%
    dplyr::rename(reads = tidyselect::all_of(colname)) %>%
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
    dplyr::mutate(chr = factor(chr, levels = (CNV_RegionsWithReads %>% pull(chr) %>% levels()))) %>%
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


#' Run HMMcopy on per-window reads for a single sample
#'
#' Given a table of per-window counts with GC and mappability, run HMMcopy
#' correction and segmentation, and return a GRanges of CNV calls for the
#' selected sample/column.
#'
#' @param CNV_RegionsWithReads A \code{data.frame} / \code{data.table} / \linkS4class{GRanges}
#'   coercible object with columns:
#'   \code{chr}, \code{start}, \code{end}, \code{width}, \code{strand}, \code{gc}, \code{map},
#'   and one sample-specific \code{reads} column (passed by name via \code{colname}).
#' @param colname Character(1). The sample column to use as read counts.
#' @param plotDir Character(1) or \code{NULL}. If provided, writes a PDF plot named \code{<colname>.pdf}.
#'
#' @return A \linkS4class{GRanges} with columns named \code{<colname>} containing
#'   per-window CNV estimates (HMM state means on the log2 scale).
#'
#' @seealso \code{\link{addHMMcopyCNV}}, \pkg{HMMcopy}
#' @family CNV
#'
#' @examples
#' # Minimal skeleton (non-running): build the required table, then call:
#' # tbl <- data.frame(
#' #   chr = rep(paste0("chr", 1:2), each = 5),
#' #   start = seq(1, by = 1e6, length.out = 10),
#' #   end   = seq(1e6, by = 1e6, length.out = 10),
#' #   width = 1e6, strand = "*",
#' #   gc = runif(10, 0.3, 0.7), map = runif(10, 0.8, 1.0),
#' #   sampleA = rpois(10, 500)
#' # )
#' # gr_cnv <- runHMMCopy(tbl, colname = "sampleA", plotDir = tempdir())
#' # gr_cnv
#' 
#' @export
plotCNVheatmap <- function(qseaSet,
                           sampleAnnotation = NULL,
                           annotationColors = NA,
                           clusterRows = TRUE){

  rowAnnot <- makeHeatmapAnnotations(qseaSet,
                                     sampleOrientation = "row",
                                     specifiedAnnotationColors = annotationColors,
                                     sampleAnnotation = {{sampleAnnotation}} ) %>%
    purrr::pluck("sample")

  chr <- qseaSet %>%
    qsea::getCNV() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = factor(seqnames, levels = gtools::mixedsort(levels(seqnames)))) %>%
    dplyr::mutate(window = paste0(seqnames, ":",start, "-",end)) %>%
    dplyr::arrange(seqnames) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("window") %>%
    dplyr::pull(seqnames)

  chr.levs <- chr %>% levels()
  chr.cols <- list(chr = rep(c("black","grey"), length(chr.levs) ) %>%
                     utils::head(length(chr.levs)) %>%
                     purrr::set_names(chr.levs))

  #Make top_annotation bar indicating chromosomes
  topAnnot <- ComplexHeatmap::HeatmapAnnotation(chr = chr, col = chr.cols, show_legend = FALSE, show_annotation_name = FALSE)

  CNVmatrix <- qseaSet %>%
    qsea::getCNV() %>%
    tibble::as_tibble(rownames = NULL) %>%
    dplyr::mutate(window = paste0(seqnames, ":",start, "-",end)) %>%
    tibble::column_to_rownames("window") %>%
    dplyr::select(tidyselect::all_of(qseaSet %>% qsea::getSampleNames())) %>%
    as.matrix() %>%
    t()

  CNVmatrix %>%
    ComplexHeatmap::Heatmap(cluster_rows = clusterRows,
                            cluster_columns = FALSE,
                            show_column_names = FALSE,
                            left_annotation = rowAnnot,
                            top_annotation = topAnnot,
                            name = "Copy number",
                            row_title = "Samples",
                            column_title = "Chromosome",
                            column_title_side = "top",
                            heatmap_legend_param = list(legend_direction = "horizontal")) %>%
    ComplexHeatmap::draw(heatmap_legend_side = "bottom",
                         annotation_legend_side = "right")

}


#' Remove (zero out) CNV data from a qseaSet
#'
#' Sets all per-window CNV values to zero, keeping genomic coordinates intact.
#' Useful for resetting CNV prior to recalculation.
#'
#' @param qseaSet A \code{qseaSet}.
#' @return A \code{qseaSet} with the CNV assay zeroed.
#' @seealso \code{\link{addHMMcopyCNV}}, \code{\link{plotCNVheatmap}}
#' @family CNV
#'
#' @examples
#' # data(exampleTumourNormal, package = "mesa")
#' # qs <- exampleTumourNormal
#' # qs <- removeCNV(qs)
#' # qsea::getCNV(qs)
#' 
#' @export
removeCNV <- function(qseaSet){

  qseaSet@cnv <- qseaSet@cnv %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(-c(seqnames, start, end, width, strand),~.*0)) %>%
    plyranges::as_granges()

  return(qseaSet)

}
