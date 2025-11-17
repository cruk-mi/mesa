#' Add per-sample CNV calls to a qseaSet using HMMcopy
#'
#' Runs a CNV pipeline on each sample: counts fragments in fixed-size genomic
#' windows, normalises for GC and mappability (HMMcopy), segments, and stores
#' per-window copy-number in the \code{qseaSet}. Supports parallel execution.
#' Only GRCh38/hg38 is currently supported.
#'
#' @param qs `qseaSet`.  
#'   Input object whose samples will receive CNV tracks.
#'
#' @param inputColumn `character(1)`.  
#'   Column in the sample table containing BAM paths (e.g. `"input_file"`).  
#'   **Default:** `"input_file"`.
#'
#' @param windowSize `integer(1)`.  
#'   CNV bin size in bp. Allowed values: `50000`, `500000`, `1000000`.  
#'   **Default:** `1000000`.
#'
#' @param fragmentLength `integer(1)` or `NULL`.  
#'   Length to extend unpaired R1 reads; if `NULL`, estimate from proper pairs.  
#'   **Default:** `NULL`.
#'
#' @param plotDir `character(1)` or `NULL`.  
#'   Directory to write per-sample PDF diagnostics; if `NULL`, no plots are saved.  
#'   **Default:** `NULL`.
#'
#' @param parallel `logical(1)`.  
#'   Use the registered **BiocParallel** backend when `TRUE`; otherwise run serially.  
#'   **Default:** `getMesaParallel()`.
#'
#' @param maxInsertSize `integer(1)`.  
#'   Maximum insert size to accept for proper pairs (bp).  
#'   **Default:** `1000`.
#'
#' @param minInsertSize `integer(1)`.  
#'   Minimum insert size to accept for proper pairs (bp).  
#'   **Default:** `50`.
#'
#' @param minReferenceLength `integer(1)`.  
#'   Minimum aligned reference span for unpaired R1 (bp).  
#'   **Default:** `30`.
#'
#' @param minMapQual `integer(1)`.  
#'   Minimum MAPQ; for pairs, either end meeting the cutoff is accepted.  
#'   **Default:** `30`.
#'
#' @param properPairsOnly `logical(1)`.  
#'   If `TRUE`, ignore unpaired R1; if `FALSE`, combine proper pairs with high-quality R1.  
#'   **Default:** `FALSE`.
#'
#' @param hmmCopyGC `GRanges` or `NULL`.  
#'   Precomputed GC content track binned at `windowSize`.  
#'   **Default:** `NULL`.
#'
#' @param hmmCopyMap `GRanges` or `NULL`.  
#'   Precomputed mappability track binned at `windowSize`.  
#'   **Default:** `NULL`.
#'
#' @details
#' Coverage is obtained (pairs + optional R1) over tiled windows, corrected with
#' HMMcopy (`HMMcopy::correctReadcount()` followed by `HMMcopy::HMMsegment()`), and
#' merged back into the `qseaSet`. You must supply GC/mappability tracks binned at
#' the same `windowSize` (or use package-provided defaults where applicable).
#' When `parallel = TRUE`, computation uses the currently registered backend
#' (see `BiocParallel::register()`).
#'
#' @return
#' The input `qseaSet` with CNV information attached, including:
#' - per-window CNV tracks accessible via downstream CNV accessors,
#' - library/sample tables updated with relevant CNV/coverage summaries,
#' - (optionally) per-sample diagnostic PDFs if `plotDir` is provided.
#'
#' @seealso
#' \code{runHMMcopy()}, 
#' [plotCNVheatmap()], 
#' \code{getBamCoveragePairedAndUnpairedR1()} (internal low-level worker),  
#' [HMMcopy::correctReadcount()], 
#' [HMMcopy::HMMsegment()], 
#' [BiocParallel::register()]
#'
#' @family CNV
#'
#' @examples
#' \dontrun{
#' data(exampleTumourNormal, package = "mesa")
#' exampleTumourNormal %>%
#'
#'     addHMMcopyCNV(
#'            inputColumn = "input_file",
#'            windowSize    = 1e6,
#'            hmmCopyGC     = gc_hg38_1000kb,   # GRanges with 'gc' column
#'            hmmCopyMap    = map_hg38_1000kb,  # GRanges with 'map' column
#'            properPairsOnly = TRUE,
#'            parallel      = FALSE,
#'            plotDir       = tempdir()
#'      )
#' 
#' }
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
    stop("Mappability track for hmmcopy must be provided.")
  }

  hmmCopyGC <- asValidGranges(hmmCopyGC)
  hmmCopyMap <- asValidGranges(hmmCopyMap)

  #TODO something with zygosity

  CNV_Regions <- qsea:::makeGenomeWindows(qsea:::getGenome(qs), as.character(qsea::getChrNames(qs)), windowSize)

  # check that the hmmcopy objects are appropriate; each CNV_Region window should
  # overlap exactly one GC/Map window. 
  overlappedRegions <- CNV_Regions %>%
    plyranges::join_overlap_left(hmmCopyGC) %>%
    plyranges::join_overlap_left(hmmCopyMap)
  
  # two conditions should be the same, but really don't want this to happen as it 
  # leads to exponentially increasing numbers of rows 
  if (any(overlappedRegions %>% duplicated()) ||
      (length(overlappedRegions) != length(CNV_Regions))) {
    
    duplicatedRegions <- utils::head(
      overlappedRegions[overlappedRegions %>% duplicated()]
    )
    
    stop(
      sprintf(
        "CNV regions overlap with multiple windows from the hmmCopyGC and/or \
hmmCopyMap objects!\nThis probably means that your window sizes do not match.\n\
Showing first affected regions:\n%s",
        print_and_capture(duplicatedRegions)
      ),
      call. = FALSE
    )
  }
  
  # this object is only used for this sanity check, not used directly again, so clean up
  rm(overlappedRegions)
  
  if (parallel) {
    BPPARAM <- BiocParallel::bpparam()
    message("Scanning up to ", BiocParallel::bpnworkers(BPPARAM), " files in parallel")
  }  else {BPPARAM <- BiocParallel::SerialParam()}

  bamOutList <- BiocParallel::bplapply(X = qsea::getSampleTable(qs) %>% dplyr::pull(inputColumn),
                                       FUN = getBamCoveragePairedAndUnpairedR1,
                                       BSgenome = qsea:::getGenome(qs),
                                       regions = CNV_Regions,
                                       fragmentLength = fragmentLength, maxInsertSize = maxInsertSize,
                                       minInsertSize = minInsertSize, minReferenceLength = minReferenceLength,
                                       minMapQual = minMapQual, properPairsOnly = properPairsOnly, BPPARAM = BPPARAM)

  names(bamOutList) <- qsea::getSampleTable(qs)$sample_name

  libraries <- bamOutList %>%
    purrr::map_dfr(~ purrr::pluck(.,"library")) %>%
    as.data.frame()

  rownames(libraries) <- qsea::getSampleTable(qs)$sample_name

  qs <- qsea:::setLibrary(qs, "input_file", libraries)

  counts <-  bamOutList %>%
    purrr::map_dfc(function(x) {
      purrr::pluck(x,"regionsGR") %>%
        tibble::as_tibble() %>%
        dplyr::pull(counts)  }
    ) %>% as.matrix()
  colnames(counts) <- qsea::getSampleTable(qs)$sample_name

  nf <- apply(X = counts, MARGIN = 2, FUN = stats::median, na.rm = TRUE)
  if (any(nf < 100, na.rm = TRUE)) {
    warning("Low coverage in CNV files. Consider larger CNV_window_size")
  }

  GenomicRanges::mcols(CNV_Regions) <- tibble::as_tibble(counts)

  CNV_RegionsWithReads <- CNV_Regions %>%
      plyranges::join_overlap_left(hmmCopyGC) %>%
      plyranges::join_overlap_left(hmmCopyMap) %>%
      tibble::as_tibble() %>%
      dplyr::rename(chr = seqnames) %>%
      dplyr::mutate(
          gc  = as.numeric(gc),
          map = as.numeric(map)
      ) %>% 
      data.table::data.table()
  

  #try saving the raw CNV data in the qseaSet.
  #TODO would be better as a slot of its own, but can't seem to manage that into the S4 class
  qs@libraries$input_counts <- CNV_RegionsWithReads

  mergedCNV <- purrr::map(qsea::getSampleTable(qs)$sample_name, ~ runHMMCopy(CNV_RegionsWithReads, ., plotDir = plotDir)) %>%
    purrr::reduce(plyranges::join_overlap_inner)

  qs <- qsea:::setCNV(qs, mergedCNV)

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
#' @param CNV_RegionsWithReads `data.frame` / `data.table` / `GRanges`.  
#'   Coercible table with columns:
#'   `chr`, `start`, `end`, `width`, `strand`, `gc`, `map`, and one
#'   sample-specific count column referenced by `colname`.
#'
#' @param colname `character(1)`.  
#'   Name of the sample count column to use (e.g., `"sampleA"`).
#'
#' @param plotDir `character(1)` or `NULL`.  
#'   Directory to write a per-sample PDF diagnostic plot named `<colname>.pdf`;  
#'   if `NULL`, no plot is produced.  
#'   **Default:** `NULL`.
#'
#' @details
#' Counts are GC/map corrected via `HMMcopy::correctReadcount()` and segmented
#' with `HMMcopy::HMMsegment()`. The result is mapped back to the input windows
#' and returned as a `GRanges`. Ensure that the binning (`width`) and covariate
#' columns (`gc`, `map`) are consistent across the table.
#'
#' @return A `GRanges` with a metadata column named `<colname>` containing
#'   per-window CNV estimates (typically HMM state means on the log2 scale).
#'
#' @seealso
#' [addHMMcopyCNV()], [HMMcopy::correctReadcount()], [HMMcopy::HMMsegment()]
#'
#' @family CNV
#'
#' @examples
#' \dontrun{
#' # Minimal synthetic example: build a table and pipe into runHMMCopy()
#' set.seed(1)
#' n <- 2000L; w <- 5e4
#' tibble::tibble(
#'     chr     = factor(rep(c("chr1","chr2"), 
#'                            each = n/2), 
#'                            levels = c("chr1","chr2")),
#'     start   = rep(seq(1, by = w, length.out = n/2), 2),
#'     end     = start + w - 1L,
#'     width   = w, strand = "*",
#'     gc      = runif(n, 0.05, 0.95),        #' wide spread => stable LOESS
#'     map     = runif(n, 0.60, 0.98),        #' avoid extremes / exact 1.0
#'     sampleA = rpois(n, 500)
#' ) %>%
#'     runHMMCopy(colname = "sampleA", plotDir = NULL)
#'}
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


#' CNV heatmap across samples
#'
#' Plot per-window copy-number estimates stored in a `qseaSet` as a heatmap
#' using **ComplexHeatmap**. Each column is a sample; rows are genome windows.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object containing per-window CNV tracks (e.g., created by
#'   [addHMMcopyCNV()] or `qsea::addCNV()`).
#'
#' @param sampleAnnotation tidyselect specification, `character()`, or `NULL`.  
#'   Sample-table columns to display as column annotations (e.g., `c("tumour","type")`
#'   or bare helpers `tumour, type`).  
#'   **Default:** `NULL`.
#'
#' @param annotationColors `list` or `NA`.  
#'   Optional mapping of levels → colours to override auto palettes, e.g.
#'   `list(tumour = c(Tumour = "firebrick4", Normal = "blue"))`.  
#'   **Default:** `NA`.
#'
#' @param clusterRows `logical(1)`.  
#'   Cluster rows (windows) before plotting.  
#'   **Default:** `TRUE`.
#'
#' @details
#' Assumes CNV values are available per window and sample. Rows are windows
#' (typically in genomic order if `clusterRows = FALSE`), columns are samples.
#' Sample annotations (and optional colours) are added when `sampleAnnotation`
#' and `annotationColors` are provided.
#'
#' @return
#' Draws a heatmap on the active device and (invisibly) returns the underlying
#' **ComplexHeatmap** object for further composition with
#' [ComplexHeatmap::draw()].
#'
#' @seealso
#' [addHMMcopyCNV()], \code{runHMMcopy()}, [ComplexHeatmap::Heatmap()],
#' [qsea::addCNV()]
#'
#' @family CNV
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Basic CNV heatmap with sample annotations
#' exampleTumourNormal %>%
#'   plotCNVheatmap(sampleAnnotation = c(tumour, type))
#'
#' # Custom annotation colours; disable row clustering
#' exampleTumourNormal %>%
#'   plotCNVheatmap(
#'     sampleAnnotation = tumour,
#'     annotationColors = list(tumour = c(Tumour = "firebrick4", Normal = "blue")),
#'     clusterRows = FALSE
#'   )
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
#' Set all per-window CNV values to zero while keeping genomic coordinates
#' and sample metadata intact. Useful for resetting CNV prior to recomputation.
#'
#' @param qseaSet `qseaSet`.  
#'   Input object whose CNV assay will be replaced by zeros.
#'
#' @return A `qseaSet` with the CNV matrix zeroed (same dimensions as before).
#'
#' @details
#' This operation preserves regions (windows) and sample metadata; only the
#' per-window CNV values are replaced with zeros. Use this to “clear” existing
#' CNV calls before running [addHMMcopyCNV()] or other CNV routines again.
#'
#' @seealso
#' [addHMMcopyCNV()], [plotCNVheatmap()], [qsea::addCNV()]
#'
#' @family CNV
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Zero out CNV and inspect the first rows
#' exampleTumourNormal %>%
#'   removeCNV() %>%
#'   getCNV() %>%
#'   utils::head()
#'
#' @export
removeCNV <- function(qseaSet){

  qseaSet@cnv <- qseaSet@cnv %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(-c(seqnames, start, end, width, strand),~.*0)) %>%
    plyranges::as_granges()

  return(qseaSet)

}
