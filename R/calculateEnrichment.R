#' Genome-wide CpG statistics (relH and GoGe)
#'
#' Compute genome-wide CpG metrics for a given BSgenome:  
#' * **relH**: relative CpG frequency (% of dinucleotides that are `"CG"`).  
#' * **GoGe**: enrichment statistic \eqn{(nCG * genome_length) / (nC * nG)}.  
#'
#' These values are used as denominators to normalise sample-level CpG enrichment.
#'
#' @param BSgenome `character(1)`  
#'   The name of a BSgenome package, e.g.
#'   `"BSgenome.Hsapiens.NCBI.GRCh38"`. The package must be installed
#'   and loadable in the current session.
#'
#' @return A `data.frame` with two numeric columns:
#' * **genome.relH** — percent of `"CG"` dinucleotides genome-wide.  
#' * **genome.GoGe** — GoGe enrichment statistic.  
#'
#' @details
#' * Chromosomes with names containing `"rand"` or `"chrUn"` are excluded.  
#' * The BSgenome indicated by `BSgenome` is loaded dynamically.  
#' * Uses \pkg{Biostrings} and \pkg{BSgenome} utilities to count CpG, C, and G
#'   occurrences across the genome.  
#'
#' @seealso
#'   [calculateCGEnrichment()],  
#'   [calculateCGEnrichmentGRanges()],  
#'   \pkg{BSgenome},  
#'   \pkg{Biostrings}
#'
#' @examples
#' # Example: compute genome-wide CpG metrics for a S. cerevisae genome. 
#' # It could be done for GRCh38 genome instead, but it will take more time to
#' # run (package must be installed)
#' if (requireNamespace("BSgenome.Scerevisiae.UCSC.sacCer3", quietly = TRUE)) {
#'   calculateGenomicCGDistribution("BSgenome.Scerevisiae.UCSC.sacCer3")
#'   }
#'   
#' @export
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
#' Compute CpG enrichment metrics (relH and GoGe) from aligned reads in a BAM
#' file. Uses \pkg{MEDIPS} to obtain fragment ranges and \pkg{Biostrings} to
#' interrogate the reference genome. Optionally exports a fragment-length
#' density plot (PDF) and a serialized RDS with the GRanges of reads.
#'
#' @param file `character(1)`  
#'   Path to the BAM file.  
#'   **Default:** `NULL` (must be supplied).
#'
#' @param BSgenome `character(1)`  
#'   Name of a BSgenome package, e.g. `"BSgenome.Hsapiens.NCBI.GRCh38"`.  
#'   For GRCh38 and hg19, precomputed distributions are cached within \pkg{mesa}
#'   for speed; otherwise, [calculateGenomicCGDistribution()] is used.  
#'   **Default:** `NULL`.
#'
#' @param exportPath `character(1)` or `NULL`  
#'   Directory in which to write a fragment-length density PDF and an RDS file
#'   containing the read GRanges. If `NULL`, no files are written.  
#'   **Default:** `NULL`.
#'
#' @param extend `integer(1)`  
#'   Passed to [MEDIPS::getGRange()]. Extension length for unpaired reads.  
#'   **Default:** `0`.
#'
#' @param shift `integer(1)`  
#'   Passed to [MEDIPS::getGRange()]. Shift applied to read positions.  
#'   **Default:** `0`.
#'
#' @param uniq `integer(1)`  
#'   Passed to [MEDIPS::getGRange()]. Minimum mapping uniqueness.  
#'   **Default:** `0`.
#'
#' @param chr.select `character()` or `NULL`  
#'   Passed to [MEDIPS::getGRange()]. Subset of chromosomes to use.  
#'   **Default:** `NULL` (all chromosomes).
#'
#' @param paired `logical(1)`  
#'   Whether BAM contains paired-end reads (passed to
#'   [MEDIPS::getPairedGRange()]).  
#'   **Default:** `TRUE`.
#'
#' @return A `data.frame` with columns:
#' * **file** — input BAM path.  
#' * **relH** — sample relH normalised by genome relH.  
#' * **GoGe** — sample GoGe normalised by genome GoGe.  
#' * **nReads** — total reads considered.  
#' * **nReadsWithoutPattern** — reads lacking `"CG"` motif.  
#' * **n100bpReads** — reads with length ≥ 100 bp.  
#' * **n100bpReadsWithoutPattern** — reads ≥ 100 bp lacking `"CG"`.  
#'
#' @details
#' Reads are scanned for `"CG"` dinucleotides. Counts are normalised against
#' genome-wide expectations (relH, GoGe).  
#' * For **GRCh38** and **hg19**, cached genomic distributions are bundled
#'   with \pkg{mesa} for efficiency.  
#' * For other genomes, [calculateGenomicCGDistribution()] is invoked.  
#'
#' @seealso
#'   [calculateCGEnrichmentGRanges()],  
#'   [calculateGenomicCGDistribution()],  
#'   \pkg{MEDIPS},  
#'   \pkg{BSgenome}
#'   
#'
#' @examples
#' if (requireNamespace("MEDIPS", quietly = TRUE) &&
#'   requireNamespace("MEDIPSData", quietly = TRUE) &&
#'   requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE) &&
#'   require("GenomicRanges", quietly = TRUE)) {
#'   
#'   calculateCGEnrichment(
#'     file       = system.file("extdata", "hESCs.Input.chr22.bam", package = "MEDIPSData"),
#'     BSgenome   = "BSgenome.Hsapiens.UCSC.hg19",  
#'     exportPath = tempdir(),
#'     paired     = FALSE                            
#'   )
#' }
#' 
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
#' # Requires MEDIPS and a BSgenome package
#' # if (requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
#' #   getCGPositions("BSgenome.Hsapiens.NCBI.GRCh38", chr.select = 22)
#' # }
getCGPositions <- function(BSgenome, chr.select){
  MEDIPS::MEDIPS.getPositions(BSgenome, "CG", chr.select)
}


#' CpG enrichment from GRanges of reads (MEDIPS-style)
#'
#' Compute CpG enrichment metrics (**relH** and **GoGe**) from a
#' [GenomicRanges::GRanges-class] of read spans, using the specified BSgenome.
#' Useful when reads are already represented as genomic ranges rather than
#' extracted directly from a BAM file.
#'
#' @param readGRanges `GRanges`  
#'   Genomic ranges representing fragments (each range = one fragment).  
#'   **Default:** `NULL` (must be supplied).
#'
#' @param BSgenome `character(1)`  
#'   Name of a BSgenome package, e.g. `"BSgenome.Hsapiens.NCBI.GRCh38"`.  
#'   The package must be installed and loadable.  
#'   **Default:** `NULL`.
#'
#' @param chr.select `character()` or `NULL`  
#'   Vector of chromosomes to restrict motif calculation.  
#'   **Default:** `NULL` (use all chromosomes).
#'
#' @return A `tibble` with one row and the following columns:
#' * **relH** — sample relH normalised by genome relH.  
#' * **GoGe** — sample GoGe normalised by genome GoGe.  
#' * **nReads** — total reads (fragments).  
#' * **nReadsWithoutPattern** — reads lacking `"CG"` motif.  
#' * **n100bpReads** — reads with length ≥ 100 bp.  
#' * **n100bpReadsWithoutPattern** — reads ≥ 100 bp lacking `"CG"`.  
#'
#' @details
#' relH and GoGe are computed by counting `"CG"` dinucleotides in the provided
#' fragments and normalising by genome-wide expectations (from
#' [calculateGenomicCGDistribution()]).  
#' Unlike [calculateCGEnrichment()], this function assumes reads are already
#' summarised as genomic ranges.
#'
#' @seealso
#'   [calculateCGEnrichment()],  
#'   [calculateGenomicCGDistribution()],  
#'   \pkg{GenomicRanges},  
#'   \pkg{BSgenome}
#'
#' @examples
#' # Runnable toy example with synthetic reads over chr1 (requires a BSgenome)
#' if (requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
#'   n <- 200
#'   gr <- GenomicRanges::GRanges(
#'     seqnames = rep(1, n),
#'     ranges   = IRanges::IRanges(
#'       start = sample(1e6:2e6, n),
#'       width = sample(80:180, n, replace = TRUE)
#'     ),
#'     strand = "*"
#'   )
#'   calculateCGEnrichmentGRanges(
#'     readGRanges = gr,
#'     BSgenome    = "BSgenome.Hsapiens.NCBI.GRCh38",
#'     chr.select  = 1
#'   )
#' }
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
#' For each sample, compute MEDIPS-style CpG enrichment metrics (**relH**, **GoGe**,
#' and read counts) from its BAM and append the results to the `qseaSet` sample
#' metadata. Supports both pulldown (default) and input libraries.
#'
#' @param qseaSet `qseaSet`
#'   Input object whose sample table contains BAM paths.
#'
#' @param exportPath `character(1)` or `NULL`
#'   Directory to export per-sample fragment-length PDFs and read-`GRanges` RDS files.
#'   If `NULL`, no files are written.  
#'   **Default:** `NULL`.
#'
#' @param nonEnrich `logical(1)`
#'   If `TRUE`, treat samples as **Input** libraries (columns appended under
#'   `@libraries$input_file`). If `FALSE`, treat as **Pulldown** (columns appended
#'   under `@libraries$file_name`).  
#'   **Default:** `FALSE`.
#'
#' @param extend `integer(1)`  
#'   Passed to [MEDIPS::getGRange()]. Extension length for unpaired reads.  
#'   **Default:** `0`.
#'
#' @param shift `integer(1)`  
#'   Passed to [MEDIPS::getGRange()]. Shift applied to read positions.  
#'   **Default:** `0`.
#'
#' @param uniq `integer(1)`  
#'   Passed to [MEDIPS::getGRange()]. Minimum mapping uniqueness.  
#'   **Default:** `0`.
#'
#' @param chr.select `character()` or `NULL`
#'   Passed to MEDIPS range extraction; subset of chromosomes to analyse.  
#'   **Default:** `NULL` (all chromosomes).
#'
#' @param paired `logical(1)`
#'   Whether BAMs are paired-end (uses [MEDIPS::getPairedGRange()]).  
#'   **Default:** `TRUE`.
#'
#' @param file_name `character(1)`
#'   Column name in the sample table holding BAM paths when `nonEnrich = FALSE`.
#'   When `nonEnrich = TRUE`, the column `input_file` is used instead.  
#'   **Default:** `"file_name"`.
#'
#' @param nCores `integer(1)`
#'   Number of parallel cores for `parallel::mclapply()`. Set to `1` for serial.  
#'   **Default:** `1`.
#'
#' @return A `qseaSet` with new MEDIPS-style metrics appended:
#'
#' * **Sample/library metrics** (added under the appropriate library frame):
#'   `relH`, `GoGe`, `nReads`, `nReadsWithoutPattern`, `n100bpReads`,
#'   `n100bpReadsWithoutPattern`.
#' * **Provenance**: any export artefacts (PDF/RDS) written to `exportPath` if provided.
#'
#' @details
#' Internally calls [calculateCGEnrichment()] per sample to derive relH/GoGe and counts,
#' then merges the results into the `@libraries` slot (`$file_name` for pulldown,
#' `$input_file` for input when `nonEnrich = TRUE`). Parallelisation uses
#' `parallel::mclapply()`; on non-Unix systems, `nCores` > 1 is ignored.
#'
#' @seealso
#'   [calculateCGEnrichment()],  
#'   [calculateCGEnrichmentGRanges()]
#'
#' @examples
#' \donttest{
#'  if (requireNamespace("MEDIPS", quietly = TRUE) &&
#'     requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
#'       
#'       bam <- system.file("extdata", "hESCs.Input.chr22.bam", package = "MEDIPSData")
#'       
#'       data.frame(
#'         sample_name = "hESC_Input_chr22",
#'         file_name   = bam,
#'         group       = "Input",
#'         input_file  = bam,
#'         stringsAsFactors = FALSE
#'       ) %>%
#'         qsea::createQseaSet("BSgenome.Hsapiens.UCSC.hg19") %>%
#'         addMedipsEnrichmentFactors(
#'           exportPath = tempdir(),
#'           chr.select = paste0("chr22"),
#'           paired     = FALSE
#'           # If your qsea returns a BSgenome object and as.character() doesn't help, uncomment:
#'           # , BSgenome_pkg = "BSgenome.Hsapiens.UCSC.hg19"
#'         ) %>%
#'         qsea::getSampleTable() %>%
#'         print()
#'     }
#' }
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
    typeString <- "non-enriched"
  } else {
    typeString <- "enriched"
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

