#' Genome-wide CpG statistics (relH and GoGe)
#'
#' Compute genome-wide CpG metrics for a given BSgenome: the relative CpG
#' frequency (\emph{relH}, percent of dinucleotides that are "CG") and the
#' GoGe statistic \eqn{(nCG * genome_length) / (nC * nG)}. These values are
#' used as denominators to normalise sample-level CpG enrichment.
#'
#' @param BSgenome Character(1). The BSgenome package name, e.g.
#'   \code{"BSgenome.Hsapiens.NCBI.GRCh38"}.
#'
#' @return A \code{data.frame} with two columns:
#' \describe{
#'   \item{\code{genome.relH}}{Percent of "CG" dinucleotides genome-wide.}
#'   \item{\code{genome.GoGe}}{GoGe enrichment statistic.}
#' }
#'
#' @details
#' Chromosomes named \code{"rand"} or \code{"chrUn"} are excluded. The function
#' loads the BSgenome indicated by \code{BSgenome} and uses \pkg{Biostrings}
#' and \pkg{BSgenome} utilities to count CG, C, and G across the genome.
#'
#' @seealso \code{\link{calculateCGEnrichment}},
#'   \code{\link{calculateCGEnrichmentGRanges}}, \pkg{BSgenome}, \pkg{Biostrings}
#'
#' @examples
#' \donttest{
#' # Example: compute genome-wide CpG metrics for GRCh38 (package must be installed)
#' if (requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
#'   calculateGenomicCGDistribution("BSgenome.Hsapiens.NCBI.GRCh38")
#' }
#' }
calculateGenomicCGDistribution <- function(BSgenome){
  dataset = eval(parse(text=paste0(BSgenome,"::", BSgenome)))
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


#' CpG enrichment from a BAM file (MEDIPS-style)
#'
#' Compute CpG enrichment metrics (relH and GoGe) from aligned reads in a BAM,
#' using \pkg{MEDIPS} to obtain fragment ranges and \pkg{Biostrings} to
#' interrogate the reference genome. Optionally exports a fragment-length
#' density plot (PDF) and a serialized RDS with the GRanges of reads.
#'
#' @param file Character(1). Path to a BAM file.
#' @param BSgenome Character(1). BSgenome package name (e.g., \code{"BSgenome.Hsapiens.NCBI.GRCh38"}).
#'   GRCh38 and hg19 distributions are cached within \pkg{mesa} for speed.
#' @param exportPath Character(1) or \code{NULL}. Directory to write a fragment-length
#'   density PDF and an RDS of the read GRanges. If \code{NULL}, no files are written.
#' @param extend,shift,uniq,chr.select,paired Arguments passed to
#'   \code{MEDIPS::getGRange()} or \code{MEDIPS::getPairedGRange()}.
#'
#' @return A \code{data.frame} with:
#' \describe{
#'   \item{\code{file}}{Input BAM path.}
#'   \item{\code{relH}}{Sample relH normalised by genome relH.}
#'   \item{\code{GoGe}}{Sample GoGe normalised by genome GoGe.}
#'   \item{\code{nReads}}{Total reads considered.}
#'   \item{\code{nReadsWithoutPattern}}{Reads lacking "CG" motif.}
#'   \item{\code{n100bpReads}}{Reads with length \eqn{\ge} 100 bp.}
#'   \item{\code{n100bpReadsWithoutPattern}}{Reads \eqn{\ge} 100 bp lacking "CG".}
#' }
#'
#' @details
#' The function counts exact "CG" matches in read sequences and normalises by
#' genome-wide CpG metrics (relH, GoGe). For supported genomes, precomputed
#' distributions bundled in \pkg{mesa} are used; otherwise,
#' \code{\link{calculateGenomicCGDistribution}} is called.
#'
#' @seealso \code{\link{calculateCGEnrichmentGRanges}},
#'   \code{\link{calculateGenomicCGDistribution}}, \pkg{MEDIPS}, \pkg{BSgenome}
#'
#' @examples
#' \donttest{
#' # Sketch of usage (requires a real BAM and installed BSgenome):
#' # calculateCGEnrichment(
#' #   file = "sample.bam",
#' #   BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38",
#' #   exportPath = tempdir(),
#' #   paired = TRUE
#' # )
#' }
#' @export
calculateCGEnrichment <- function(file = NULL, BSgenome = NULL, exportPath = NULL,
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
  } else if (BSgenome == "BSgenome.Hsapiens.UCSC.hg19") {
    genomicDistribution <- mesa::BSgenome.Hsapiens.UCSC.hg19.CpG.distribution
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


#' Genomic positions of a motif (CG) from MEDIPS
#'
#' Convenience wrapper that returns genomic positions of the \code{"CG"} motif
#' for the specified BSgenome and chromosomes, via \code{MEDIPS::MEDIPS.getPositions()}.
#'
#' @param BSgenome Character(1). BSgenome package name.
#' @param chr.select Character vector of chromosome names to include (e.g., \code{paste0("chr", 1:22)}).
#'
#' @return A \linkS4class{GRanges} of motif positions.
#'
#' @seealso \code{\link{calculateCGEnrichment}}, \code{\link{calculateCGEnrichmentGRanges}}, \pkg{MEDIPS}
#'
#' @examples
#' \donttest{
#' # Requires MEDIPS and a BSgenome package
#' # if (requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
#' #   gr <- getCGPositions("BSgenome.Hsapiens.NCBI.GRCh38", chr.select = paste0("chr", 1:2))
#' #   gr
#' # }
#' }
getCGPositions <- function(BSgenome, chr.select){
  MEDIPS::MEDIPS.getPositions(BSgenome, "CG", chr.select)
}


#' CpG enrichment from GRanges of reads (MEDIPS-style)
#'
#' Compute relH and GoGe CpG enrichment metrics from a \linkS4class{GRanges} of
#' read spans, using the specified BSgenome. Useful when reads are already
#' represented as genomic ranges rather than a BAM file.
#'
#' @param readGRanges \linkS4class{GRanges}. Each range represents a fragment.
#' @param BSgenome Character(1). BSgenome package name (e.g., \code{"BSgenome.Hsapiens.NCBI.GRCh38"}).
#' @param chr.select Character vector of chromosomes to restrict the motif positions calculation.
#'
#' @return A \code{tibble} with columns:
#' \code{relH}, \code{GoGe}, \code{nReads}, \code{nReadsWithoutPattern},
#' \code{n100bpReads}, \code{n100bpReadsWithoutPattern}.
#'
#' @seealso \code{\link{calculateCGEnrichment}}, \code{\link{calculateGenomicCGDistribution}}
#'
#' @examples
#' \donttest{
#' # Runnable toy example with synthetic reads over chr1 (requires a BSgenome)
#' if (requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 200
#'   gr <- GenomicRanges::GRanges(
#'     seqnames = rep("chr1", n),
#'     ranges   = IRanges::IRanges(
#'       start = sample(1e6:2e6, n),
#'       width = sample(80:180, n, replace = TRUE)
#'     ),
#'     strand = "*"
#'   )
#'   calculateCGEnrichmentGRanges(
#'     readGRanges = gr,
#'     BSgenome    = "BSgenome.Hsapiens.NCBI.GRCh38",
#'     chr.select  = "chr1"
#'   )
#' }
#' }
#' 
#' @export
calculateCGEnrichmentGRanges <- function(readGRanges = NULL, BSgenome = NULL, chr.select = NULL){

  if (!requireNamespace("MEDIPS", quietly = TRUE)) {
    stop(
      "Package \"MEDIPS\" must be installed to use this function.",
      call. = FALSE
    )
  }

  dataset = eval(parse(text=paste0(BSgenome,"::", BSgenome)))

  chromosomes = gtools::mixedsort(unique(GenomeInfoDb::seqlevels(readGRanges)))

  chr_lengths = as.numeric(GenomeInfoDb::seqlengths(dataset)[chromosomes])

  if (all(is.na(GenomeInfoDb::seqlengths(readGRanges)))) {
    suppressWarnings(GenomeInfoDb::seqinfo(readGRanges) <- GenomeInfoDb::seqinfo(BSgenome::getBSgenome(BSgenome))[GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(readGRanges))])
    readGRanges <- IRanges::trim(readGRanges)
  }

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
  } else  if (BSgenome == "BSgenome.Mmusculus.UCSC.mm10") {
    genomicDistribution <- mesa::BSgenome.Mmusculus.UCSC.mm10.CpG.distribution
  } else if (BSgenome == "BSgenome.Hsapiens.UCSC.hg19") {
    genomicDistribution <- mesa::BSgenome.Hsapiens.UCSC.hg19.CpG.distribution
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


#' Add MEDIPS-style enrichment metrics to a qseaSet sample table
#'
#' For each sample, compute MEDIPS-style CpG enrichment metrics (relH, GoGe,
#' and counts) from its BAM, and append the results to the sample table within
#' the \code{qseaSet}. Supports both pulldown (default) and input libraries.
#'
#' @param qseaSet A \code{qseaSet}.
#' @param exportPath Character(1) or \code{NULL}. Directory to export per-sample
#'   fragment-length PDFs and read GRanges RDS files (optional).
#' @param nonEnrich Logical(1). If \code{TRUE}, treat samples as Input; else Pulldown.
#'   Determines which columns are appended.
#' @param extend,shift,uniq,chr.select,paired Passed to \pkg{MEDIPS} range extraction.
#' @param file_name Character(1). Column name in the sample table holding BAM paths
#'   (used when \code{nonEnrich = FALSE}). For Input, the column \code{input_file} is used.
#' @param nCores Integer(1). Number of parallel cores for \code{parallel::mclapply()}.
#'
#' @return The input \code{qseaSet}, with new columns appended to the corresponding
#'   \code{@libraries} slot (either \code{$file_name} or \code{$input_file} frame),
#'   including relH/GoGe and read counts.
#'
#' @seealso \code{\link{calculateCGEnrichment}}, \code{\link{calculateCGEnrichmentGRanges}}
#'
#' @examples
#' \donttest{
#' # Sketch (requires BAMs and a BSgenome):
#' # data(exampleTumourNormal, package = "mesa")
#' # qs <- exampleTumourNormal
#' # qs <- addMedipsEnrichmentFactors(
#' #   qs, exportPath = tempdir(), nonEnrich = FALSE,
#' #   file_name = "file_name", nCores = 1,
#' #   paired = TRUE, chr.select = paste0("chr", 1:22)
#' # )
#' # head(qsea::getSampleTable(qs))
#' }
#' 
#' @export
addMedipsEnrichmentFactors <- function(qseaSet, exportPath = NULL, nonEnrich = FALSE,
                                       extend = 0, shift = 0, uniq = 0,
                                       chr.select = NULL, paired = TRUE,
                                       file_name = "file_name",
                                       nCores = 1){

  BSgenome = qseaSet %>% qsea:::getGenome()

  if (!requireNamespace("MEDIPS", quietly = TRUE)) {
    stop(
      "Package \"MEDIPS\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if(nonEnrich){
    typeString <- "Pulldown"
  } else {
    typeString <- "Input"
  }

  message(glue::glue("Adding Medips Enrichment factors to {length(getSampleNames(qseaSet))} {typeString} samples, using {nCores} cores."))

  if(!nonEnrich){

  colsToCheck <- c("relH","GoGe","nReads","nReadsWithoutPattern","n100bpReadsWithoutPattern")

  if (any(colsToCheck %in% colnames(qsea::getSampleTable(qseaSet)))) {
    stop(glue::glue("Column {colsToCheck[colsToCheck %in% colnames(qsea::getSampleTable(qseaSet))]} already in sampleTable!
                    "))
     }
  } else{

    colsToCheck <- c("input_relH","input_GoGe","input_nReads","input_nReadsWithoutPattern","input_n100bpReadsWithoutPattern")

    if (any(colsToCheck %in% colnames(qsea::getSampleTable(qseaSet)))) {
      stop(glue::glue("Column {colsToCheck[colsToCheck %in% colnames(qsea::getSampleTable(qseaSet))]} already in sampleTable!
                    "))
    }

  }

  if(!nonEnrich){
    fileNames <- qsea::getSampleTable(qseaSet) %>% dplyr::pull(file_name)
  } else {
    fileNames <- qsea::getSampleTable(qseaSet) %>% dplyr::pull(input_file)
  }


  enrichData <- parallel::mclapply(fileNames,
                                        function(x){
                                          calculateCGEnrichment(x, BSgenome = BSgenome, exportPath = exportPath,
                                                                 extend = extend, shift = shift, uniq = uniq,
                                                                 chr.select = chr.select, paired = paired)
                                        }, mc.cores = nCores) %>%
    do.call(rbind, . )

  if(!nonEnrich){
  qseaSet@libraries$file_name <- qseaSet@libraries$file_name %>%
    cbind(dplyr::select(enrichData,-file))

  } else{
    qseaSet@libraries$input_file <- qseaSet@libraries$input_file %>%
      cbind(dplyr::select(enrichData,-file))
  }

  return(qseaSet)
}

