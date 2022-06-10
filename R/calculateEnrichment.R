#' This function takes a BSGenome and returns a data frame with the relH and GoGe parameters used for calculating enrichment
#'
#' @param BSgenome The BSgenome object.
#' @return A data frame with two entries, genome.relH and genome.GoGe
#' @export
calculateGenomicCGDistribution <- function(BSgenome){
  dataset = get(ls(paste("package:", BSgenome, sep = "")))
  CG <- Biostrings::DNAStringSet("CG")
  pdict0 <- Biostrings::PDict(CG)
  params <- methods::new("BSParams", X = dataset, FUN = Biostrings::countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
  genome.CG = sum(BSgenome::bsapply(params, pdict = pdict0))
  params <- methods::new("BSParams", X = dataset, FUN = Biostrings::alphabetFrequency, simplify = TRUE, exclude = c("rand", "chrUn"))
  alphabet = BSgenome::bsapply(params)
  genome.l = sum(as.numeric(alphabet))
  genome.C = as.numeric(sum(alphabet[2,]))
  genome.G = as.numeric(sum(alphabet[3,]))
  genome.relH = genome.CG/genome.l * 100
  genome.GoGe = (genome.CG * genome.l)/(genome.C * genome.G)
  return(data.frame(genome.relH = genome.relH, genome.GoGe = genome.GoGe))
}


#' This function takes a path to a bam file and calculates the enrichment scores.
#'
#' @param file The path to the bam file
#' @param BSgenome Which BSGenome to use (GRCh38 is cached for speed)
#' @param exportPath Location to save a plot and data object into.
#' @param extend,shift,uniq,chr.select,paired Options to feed into medips::getGRange or medips::getPairedGRange
#' @return A data frame containing the relH, GoGe and number of reads values for the samples
#' @export
calculateCpGEnrichment <- function(file = NULL, BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38", exportPath = NULL,
                                   extend = 0, shift = 0, uniq = 0,
                                   chr.select = NULL, paired = TRUE){

  if (!requireNamespace("MEDIPS", quietly = TRUE)) {
    stop(
      "Package \"MEDIPS\" must be installed to use this function.",
      call. = FALSE
    )
  }

  dataset = eval(parse(text=paste0(BSgenome,"::", BSgenome)))

  ## Read region file
  fileName = basename(file)
  path = dirname(file)
  if (path == "") {path = getwd()}
  if (!fileName %in% dir(path)) {stop(paste0("File", fileName, " not found in", path))}

  if (!paired) {GRange.Reads = MEDIPS::getGRange(fileName, path, extend, shift, chr.select, dataset, uniq, simpleCigar = FALSE)} else
  {GRange.Reads = MEDIPS::getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq, simpleCigar = FALSE)}

  ## Sort chromosomes
  if (length(unique(GenomeInfoDb::seqlevels(GRange.Reads))) > 1) {chromosomes = gtools::mixedsort(unique(GenomeInfoDb::seqlevels(GRange.Reads)))}
  if (length(unique(GenomeInfoDb::seqlevels(GRange.Reads))) == 1) {chromosomes = unique(GenomeInfoDb::seqlevels(GRange.Reads))}

  chr_lengths = as.numeric(GenomeInfoDb::seqlengths(dataset)[chromosomes])

  IRanges::ranges(GRange.Reads) <- IRanges::restrict(IRanges::ranges(GRange.Reads), +1)

  ##Calculate CpG density for regions
  total = length(chromosomes)

  readsChars <- unlist(Biostrings::getSeq(dataset, GRange.Reads, as.character = TRUE))

  # Faster to use stringr, as we are looking for exact matches.
  regions.CG = sum(stringr::str_count(readsChars, stringr::fixed("CG")))
  regions.C  = sum(stringr::str_count(readsChars, stringr::fixed("C")))
  regions.G  = sum(stringr::str_count(readsChars, stringr::fixed("G")))
  all.genomic = sum(stringr::str_length(readsChars))

  nReads <- length(readsChars)

  regions.relH = as.numeric(regions.CG) / as.numeric(all.genomic) * 100
  regions.GoGe = (as.numeric(regions.CG) * as.numeric(all.genomic)) / (as.numeric(regions.C) * as.numeric(regions.G))

  if (BSgenome == "BSgenome.Hsapiens.NCBI.GRCh38") {
    genomicDistribution <- mesa::BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution
  } else {
    genomicDistribution <- calculateGenomicCGDistribution(BSgenome)
  }

  genome.relH <- genomicDistribution$genome.relH
  genome.GoGe <- genomicDistribution$genome.GoGe

  enrichment.score.relH = regions.relH/genome.relH
  enrichment.score.GoGe = regions.GoGe/genome.GoGe

  if (!is.null(exportPath)) {
    GRange.Reads %>%
      BiocGenerics::width() %>%
      dplyr::as_tibble() %>%
      ggplot2::ggplot(ggplot2::aes(x = value)) +
      ggplot2::geom_density(color = "black") +
      ggplot2::theme_bw() +
      ggplot2::labs(xlab = "Fragment Length",
                    ylab = "Density",
                    title = "Fragment Length Distribution",
                    subtitle = stringr::str_remove(fileName, ".bam"))

    ggplot2::ggsave(file = stringr::str_replace(file.path(exportPath, fileName), ".bam", ".pdf"))

    saveRDS(GRange.Reads, file = stringr::str_replace(file.path(exportPath, fileName), ".bam", ".rds"))
  }

  # if (BSgenome == "BSgenome.Hsapiens.NCBI.GRCh38") {
  #   genomeCGranges <- CGPatternGRCh38
  # } else {
  #   genomeCGranges <- MEDIPS::MEDIPS.getPositions("BSgenome.Hsapiens.NCBI.GRCh38", "CG", chr.select)
  # }
  genomeCGranges <- getCGPositions(BSgenome, chr.select)


  numWithoutPattern <- GRange.Reads %>%
    plyranges::filter_by_non_overlaps(genomeCGranges) %>%
    length()

  numReads100bp <- GRange.Reads %>%
    dplyr::filter(width >= 100) %>%
    length()

  numWithoutPatternOver100bp <- GRange.Reads %>%
    dplyr::filter(width >= 100) %>%
    plyranges::filter_by_non_overlaps(genomeCGranges) %>%
    length()

  gc()
  return(data.frame(file = file,
                    relH = enrichment.score.relH,
                    GoGe = enrichment.score.GoGe,
                    nReads = nReads,
                    nReadsWithoutPattern = numWithoutPattern,
                    n100bpReads = numReads100bp,
                    n100bpReadsWithoutPattern = numWithoutPatternOver100bp
  ))
}
#' This function takes a GRanges object and calculates the enrichment scores.
#'
#' @param BSgenome Which BSgenome to use
#' @param chr.select Which chromosomes to use
#' @return A data frame containing the relH, GoGe and number of reads values for the samples
#' @export
getCGPositions <- function(BSgenome, chr.select){
  MEDIPS::MEDIPS.getPositions("BSgenome.Hsapiens.NCBI.GRCh38", "CG", chr.select)
}

#' This function takes a GRanges object and calculates the enrichment scores.
#'
#' @param readGRanges A GRanges object containing the span of each fragment of DNA
#' @param BSgenome Which BSgenome to use (GRCh38 is cached for speed)
#' @param chr.select Which chromosomes to use in global calculation
#' @return A data frame containing the relH, GoGe and number of reads values for the samples
#' @export
calculateCpGEnrichmentGRanges <- function(readGRanges = NULL, BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38", chr.select = NULL){

  if (!requireNamespace("MEDIPS", quietly = TRUE)) {
    stop(
      "Package \"MEDIPS\" must be installed to use this function.",
      call. = FALSE
    )
  }

  dataset = eval(parse(text=paste0(BSgenome,"::", BSgenome)))

  chromosomes = gtools::mixedsort(unique(GenomeInfoDb::seqlevels(readGRanges)))

  chr_lengths = as.numeric(GenomeInfoDb::seqlengths(dataset)[chromosomes])

  IRanges::ranges(readGRanges) <- IRanges::restrict(IRanges::ranges(readGRanges), +1)

  ##Calculate CpG density for regions
  total = length(chromosomes)

  readsChars <- unlist(Biostrings::getSeq(dataset, readGRanges, as.character = TRUE))

  # Faster to use stringr, as we are looking for exact matches.
  regions.CG = sum(stringr::str_count(readsChars, stringr::fixed("CG")))
  regions.C  = sum(stringr::str_count(readsChars, stringr::fixed("C")))
  regions.G  = sum(stringr::str_count(readsChars, stringr::fixed("G")))
  all.genomic = sum(stringr::str_length(readsChars))

  regions.relH = as.numeric(regions.CG) / as.numeric(all.genomic) * 100
  regions.GoGe = (as.numeric(regions.CG) * as.numeric(all.genomic)) / (as.numeric(regions.C) * as.numeric(regions.G))

  if (BSgenome == "BSgenome.Hsapiens.NCBI.GRCh38") {
    genomicDistribution <- mesa::BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution
  } else {
    genomicDistribution <- calculateGenomicCGDistribution(BSgenome)
  }

  genome.relH <- genomicDistribution$genome.relH
  genome.GoGe <- genomicDistribution$genome.GoGe

  enrichment.score.relH = regions.relH/genome.relH
  enrichment.score.GoGe = regions.GoGe/genome.GoGe

  genomeCGranges <- getCGPositions(BSgenome, chr.select)

  numReads <- length(readsChars)

  numWithoutPattern <- readGRanges %>%
    plyranges::filter_by_non_overlaps(genomeCGranges) %>%
    length()

  numReads100bp <- readGRanges %>%
    dplyr::filter(width >= 100) %>%
    length()

  numWithoutPatternOver100bp <- readGRanges %>%
    dplyr::filter(width >= 100) %>%
    plyranges::filter_by_non_overlaps(genomeCGranges) %>%
    length()

  gc()

  return(tibble::tibble(relH = enrichment.score.relH,
                        GoGe = enrichment.score.GoGe,
                        nReads = numReads,
                        nReadsWithoutPattern = numWithoutPattern,
                        n100bpReads = numReads100bp,
                        n100bpReadsWithoutPattern = numWithoutPatternOver100bp)
  )
}

#' This function takes a qseaSet and adds the Medips-style enrichment scores to the object sampleTable.
#' @param qseaSet The cutoff to use on the windows for each sample
#' @param BSgenome Name of BSgenome package
#' @param exportPath Folder to export files
#' @param extend,shift,uniq,chr.select,paired Options to feed into medips::getGRange or medips::getPairedGRange
#' @param file_name Column of the sampleTable which contains the file_name.
#' @param nCores Number of cores to use for parallelisation
#' @return A qseaSet supplemented with more columns on the sampleTable, related to enrichment.
#' @export
addMedipsEnrichmentFactors <- function(qseaSet, BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38", exportPath = NULL,
                                       extend = 0, shift = 0, uniq = 0,
                                       chr.select = NULL, paired = TRUE,
                                       file_name = "file_name",
                                       nCores = 1){

  if (!requireNamespace("MEDIPS", quietly = TRUE)) {
    stop(
      "Package \"MEDIPS\" must be installed to use this function.",
      call. = FALSE
    )
  }

  message(glue::glue("Adding Medips Enrichment factors to {length(getSampleNames(qseaSet))} samples, using {nCores} cores."))

  colsToCheck <- c("relH","GoGe","nReads","nReadsWithoutPattern","n100bpReadsWithoutPattern")

  if (any(colsToCheck %in% colnames(qsea::getSampleTable(qseaSet)))) {
    stop(glue::glue("Column {colsToCheck[colsToCheck %in% colnames(qsea::getSampleTable(qseaSet))])} already in sampleTable!"))
  }

  pulldownFileNames <- qsea::getSampleTable(qseaSet) %>% dplyr::pull(file_name)

  enrichDataMeCap <- parallel::mclapply(pulldownFileNames,
                                        function(x){
                                          calculateCpGEnrichment(x, BSgenome, exportPath = exportPath,
                                                                 extend = extend, shift = shift, uniq = uniq,
                                                                 chr.select = chr.select, paired = paired)
                                        }, mc.cores = nCores) %>%
    do.call(rbind, . )

  # inputFileNames <- getSampleTable(qseaSet) %>% pull(input_file)
  # enrichDataInput <- parallel::mclapply(inputFileNames ,
  #                                  function(x){
  #                   calculateCpGEnrichment(x, BSgenome, exportPath = exportPath,
  #                                      extend = extend, shift = shift, uniq = uniq,
  #                                      chr.select = chr.select, paired = paired)
  #                                    }, mc.cores = nCores) %>%
  #   do.call(rbind, . )

  qseaSet@sampleTable <- qsea::getSampleTable(qseaSet) %>%
    dplyr::bind_cols(dplyr::select(enrichDataMeCap,-file))

  return(qseaSet)
}

