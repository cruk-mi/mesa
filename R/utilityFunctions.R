utils::globalVariables(c("chr", "seqnames", "start","end",".","annotation", "value","library_factor",
                         "ROI_start","ROI_end","value","window_start","window_end","window","sample_name",
                         "width","group","relH","CpG_density","strand","state","gene_name","map","CNV",
                         "EMSEMBL","SYMBOL","GENENAME","geneChr","geneStart","geneEnd","geneLength","n",
                         "adjPval","group1","group2","sample1","sample2",".up","hyperStableFractionp8","hyperStableFractionp9","hyperStableEdgar",
                         "isProperPair","isUnmappedQuery","isSupplementaryAlignment","isDuplicate","hasUnmappedMate",
                       "isNotPassingQualityControls","rname","pos","isize","MQ","mapq","isFirstMateRead","isPaired",
                       "cigar", ".rowID", "feature","annoShort","type",".comparison",".ext",".value","total_fragments",
                       "isSecondaryAlignment","ROI_ID","ID","ENSEMBL","deltaBeta",
                       "map_hg38_1000kb", "map_hg38_10kb", "map_hg38_50kb", "map_hg38_500kb",
                       "gc_hg38_1000kb","gc_hg38_50kb","gc_hg38_500kb","gc_hg38_10kb","score",
                       "log2FC","group1new","sample2new","counts","tumour","rowIndex","nUp","nDown","landscape","shortAnno","nOverCutoff",
                       "afterOverBackNum","initialOverBackNum", "qname","inOut","nSign","nStrands",
                       "chromosome_name","start_position","end_position", "input_file",
                       "topVarNumInput", "windowSdName", "pcaName",
                       "fnValue", "group2new","gap"))

#' This function takes a genomic ranges object or data frame, and lifts over from hg19 to hg38.
#'
#' @param grOrDf Either a genomic ranges object or a data frame (with seqnames column), with or without chr in the seqnames
#' @return A genomic ranges object object with coordinates in hg38.
#' @export

liftOverHg19 <- function(grOrDf){

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "Package \"rtracklayer\" must be installed to use this function.",
      call. = FALSE
    )
  }

  regions <- grOrDf %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = ifelse(stringr::str_detect(seqnames,"chr"), seqnames, paste0("chr", seqnames))) %>%
    plyranges::as_granges()

  message("Performing coordinate liftover from hg19 to hg38")
  liftover_result <- rtracklayer::liftOver(regions, mesa::hg19ToHg38.over.chain)
  GenomeInfoDb::genome(liftover_result) <- "hg38"

  return(unlist(liftover_result))
}

.getMesaGenome <- function() {
  getOption("mesa_genome", default="hg38")
}

#' This function sets a default genome for the annotateWindows function to use
#'
#' Currently supporting "hg38" and "GRCh38", this sets the TxDb and annoDb options to appropriate files.
#'
#' @param genome String corresponding to a genome.
#' @return None
#' @export

setMesaGenome <- function(genome){
  options("mesa_genome" = genome)

  return(invisible(TRUE))
}

#' Toggle whether to use parallelisation or not inside various functions in mesa
#'
#' @param nCores How many cores to set. Assumes you want to use MulticoreParam, to use another architecture use the useParallel argument and specify manually.
#' @param useParallel Boolean denoting whether or not to use parallelisation for various functions in mesa
#' @param verbose Boolean to determine whether to print messages or not
#' @return None
#' @export

setMesaParallel <- function(nCores = NULL, useParallel = FALSE, verbose = TRUE){

  if(!is.null(nCores) && nCores > 1) {
    useParallel = TRUE
    BiocParallel::register(BiocParallel::MulticoreParam(workers = nCores))
  }

  options("mesa_parallel" = useParallel)

  if(useParallel & verbose){
    message(glue::glue("Parallelisation turned on for the functions in the mesa package, currently using {BiocParallel::bpworkers()} cores.
                       Control the number of cores by calling BiocParallel::register."))
  } else if (verbose) {
    message("Parallelisation turned off for all functions in the mesa package.")
  }

  return(invisible(useParallel))
}

getMesaParallel <- function(verbose = FALSE) {

  if(verbose) {
    if (is.null(getOption("mesa_parallel"))) {
      message("Parallelisation not set, using serial evaluation. Call setMesaParallel to set for all mesa functions.")
    } else if (getOption("mesa_parallel")) {
      message(glue::glue("Using parallelisation, over {BiocParallel::bpworkers()} cores."))
    } else {
      message("Using serial evaluation")
    }
  }

  if (is.null(getOption("mesa_parallel"))) {
    return(FALSE)
  } else {
    return(getOption("mesa_parallel"))
  }

}

expect_no_error <- function(object) {
  testthat::expect_error({{ object }}, NA)
}
