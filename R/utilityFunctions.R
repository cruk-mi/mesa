utils::globalVariables(c("chr", "seqnames", "start","end",".","annotation", "value","library_factor",
                         "ROI_start","ROI_end","value","window_start","window_end","window","sample_name",
                         "width","group","relH","CpG_density","strand","state","gene_name","map","CNV",
                         "EMSEMBL","SYMBOL","GENENAME","geneChr","geneStart","geneEnd","geneLength","n",
                         "adjPval","sample1","sample2",".up","hyperStableFractionp8","hyperStableFractionp9","hyperStableEdgar",
                         "isProperPair","isUnmappedQuery","isSupplementaryAlignment","isDuplicate","hasUnmappedMate",
                       "isNotPassingQualityControls","rname","pos","isize","MQ","mapq","isFirstMateRead","isPaired",
                       "cigar", ".rowID", "feature","annoShort","type",".comparison",".ext",".value","total_fragments",
                       "isSecondaryAlignment","ROI_ID","ID","ENSEMBL","betaDelta",
                       "map_hg38_1000kb", "map_hg38_10kb", "map_hg38_50kb", "map_hg38_500kb",
                       "gc_hg38_1000kb","gc_hg38_50kb","gc_hg38_500kb","gc_hg38_10kb","score",
                       "log2FC","sample1new","sample2new","counts","tumour","rowIndex","nUp","nDown","landscape","shortAnno","nOverCutoff",
                       "afterOverBackNum","initialOverBackNum", "qname","inOut","nSign","nStrands",
                       "chromosome_name","start_position","end_position"))

#' This function takes a genomic ranges object or data frame, and lifts over from hg19 to hg38.
#'
#' @param grOrDf Either a genomic ranges object or a data frame (with seqnames column), with or without chr in the seqnames
#' @return A genomic ranges object object with coordinates in hg38.
#' @export

liftOverhg19 <- function(grOrDf){

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "Package \"rtracklayer\" must be installed to use this function.",
      call. = FALSE
    )
  }

  regions <- grOrDf %>%
    tibble::as_data_frame() %>%
    dplyr::mutate(seqnames = ifelse(stringr::str_detect(seqnames,"chr"), seqnames, paste0("chr", seqnames))) %>%
    plyranges::as_granges()

  message("Performing coordinate liftover from hg19 to hg38")
  liftover_result <- rtracklayer::liftOver(regions, mesa::hg19ToHg38.over.chain)
  GenomeInfoDb::genome(liftover_result) <- "hg38"

  return(unlist(liftover_result))
}

.getGenome <- function() {
  getOption("mesa_genome", default="hg38")
}


expect_no_error <- function(object) {
  testthat::expect_error({{ object }}, NA)
}