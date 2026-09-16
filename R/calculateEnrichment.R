#' Genome-wide CpG statistics (relH and GoGe)
#'
#' Compute genome-wide CpG metrics for a given BSgenome:
#' * **relH**: relative CpG frequency (% of dinucleotides that are `"CG"`).
#' * **GoGe**: enrichment statistic \eqn{(nCG * genome_length) / (nC * nG)}.
#'
#' These values are used as denominators to normalise sample-level CpG
#' enrichment.
#'
#' @param BSgenome `character(1)`
#' The name of a BSgenome package, e.g. `"BSgenome.Hsapiens.NCBI.GRCh38"`. The
#' package must be installed and loadable in the current session.
#'
#' @return A `data.frame` with two numeric columns:
#' * **genome.relH** — percent of `"CG"` dinucleotides genome-wide.
#' * **genome.GoGe** — GoGe enrichment statistic.
#'
#' @details
#' * Chromosomes with names containing `"rand"` or `"chrUn"` are excluded.
#' * The BSgenome indicated by `BSgenome` is loaded dynamically.
#' * Uses \pkg{Biostrings} and \pkg{BSgenome} utilities to count CpG, C, and G
#' occurrences across the genome.
#'
#' @seealso
#' [calculateCGEnrichment()], [calculateCGEnrichmentGRanges()], \pkg{BSgenome},
#' \pkg{Biostrings}
#'
#' @examples
#' # Example: compute genome-wide CpG metrics for a S. cerevisae genome.
#' # It could be done for GRCh38 genome instead, but it will take more time to
#' # run (package must be installed)
#' if (requireNamespace("BSgenome.Scerevisiae.UCSC.sacCer3", quietly = TRUE)) {
#'     calculateGenomicCGDistribution("BSgenome.Scerevisiae.UCSC.sacCer3")
#' }
#'
#' @export
calculateGenomicCGDistribution <- function(BSgenome) {
    dataset <- eval(parse(text = paste0(BSgenome, "::", BSgenome)))
    CG <- Biostrings::DNAStringSet("CG")
    pdict0 <- Biostrings::PDict(CG)
    params <- methods::new(
        "BSParams", X = dataset, FUN = Biostrings::countPDict,
        simplify = TRUE, exclude = c("rand", "chrUn")
    )
    genome.CG <- sum(BSgenome::bsapply(params, pdict = pdict0))
    params <- methods::new(
        "BSParams", X = dataset, FUN = Biostrings::alphabetFrequency,
        simplify = TRUE, exclude = c("rand", "chrUn")
    )
    alphabet <- BSgenome::bsapply(params)
    genome.l <- sum(as.numeric(alphabet))
    genome.C <- as.numeric(sum(alphabet[2, ]))
    genome.G <- as.numeric(sum(alphabet[3, ]))
    genome.relH <- genome.CG / genome.l * 100
    genome.GoGe <- (genome.CG * genome.l) / (genome.C * genome.G)
    return(data.frame(genome.relH = genome.relH, genome.GoGe = genome.GoGe))
}


#' CpG enrichment from a BAM file (MEDIPS-style)
#'
#' Compute CpG enrichment metrics (relH and GoGe) from aligned reads in a BAM
#' file. Reads are imported directly with \pkg{Rsamtools}, and
#' \pkg{Biostrings} is used to interrogate the reference genome. Optionally
#' exports a fragment-length density plot (PDF) and a serialized RDS with the
#' GRanges of reads.
#'
#' @param file `character(1)`
#' Path to the BAM file.
#'   **Default:** `NULL` (must be supplied).
#'
#' @param BSgenome `character(1)`
#' Name of a BSgenome package, e.g. `"BSgenome.Hsapiens.NCBI.GRCh38"`. For
#' GRCh38 and hg19, precomputed distributions are cached within \pkg{mesa} for
#' speed; otherwise, [calculateGenomicCGDistribution()] is used.
#'   **Default:** `NULL`.
#'
#' @param exportPath `character(1)` or `NULL`
#' Directory in which to write a fragment-length density PDF and an RDS file
#' containing the read GRanges. If `NULL`, no files are written.
#'   **Default:** `NULL`.
#'
#' @param extend `integer(1)`
#' Extension length for single-end reads (used only when `paired = FALSE`):
#' reads shorter than this are lengthened to it in the 5'->3' direction, and
#' reads already longer are left unchanged. Unused for paired reads.
#'   **Default:** `0`.
#'
#' @param shift `integer(1)`
#' Strand-aware offset applied to read positions, used only when
#' `paired = FALSE`. Unused for paired reads, where the true fragment
#' position is known.
#'   **Default:** `0`.
#'
#' @param uniq `numeric(1)`
#' How to handle duplicate fragments, following \pkg{MEDIPS}:
#' * `0` — keep every read.
#' * `1` — keep at most one read per genomic location and strand.
#' * a p-value in `(0, 1)` — keep at most `qpois(1 - uniq, rate)` reads per
#' location, where `rate` is the mean read depth per base.
#'
#' Any other value, and any logical, is an error.
#'   **Default:** `0`.
#'
#' @param chr.select `character()` or `NULL`
#' Subset of chromosomes to use. The BAM index is used to restrict the scan
#' when one is present; otherwise the whole file is scanned and the selection
#' applied afterwards.
#'   **Default:** `NULL` (all chromosomes).
#'
#' @param paired `logical(1)`
#' Whether the BAM contains paired-end reads. When `TRUE`, properly paired
#' fragments are reconstructed from the template length; when `FALSE`, single
#' reads are used (optionally extended via `extend`).
#'   **Default:** `TRUE`.
#'
#' @return A `data.frame` with columns:
#' * **file** — input BAM path.
#' * **relH** — sample relH normalised by genome relH.
#' * **GoGe** — sample GoGe normalised by genome GoGe.
#' * **nReads** — total reads considered.
#' * **nReadsWithoutPattern** — reads lacking `"CG"` motif.
#' * **n100bpReads** — reads with length >= 100 bp.
#' * **n100bpReadsWithoutPattern** — reads >= 100 bp lacking `"CG"`.
#'
#' @details
#' Reads are scanned for `"CG"` dinucleotides. Counts are normalised against
#' genome-wide expectations (relH, GoGe).
#' * For **GRCh38** and **hg19**, cached genomic distributions are bundled
#' with \pkg{mesa} for efficiency.
#' * For other genomes, [calculateGenomicCGDistribution()] is invoked.
#'
#' @seealso
#' [calculateCGEnrichmentGRanges()], [calculateGenomicCGDistribution()],
#' \pkg{BSgenome}
#'
#'
#' @examples
#' if (requireNamespace("MEDIPSData", quietly = TRUE) &&
#'     requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
#'     calculateCGEnrichment(
#'         file = system.file(
#'             "extdata",
#'             "hESCs.Input.chr22.bam",
#'             package = "MEDIPSData"
#'         ),
#'         BSgenome   = "BSgenome.Hsapiens.UCSC.hg19",
#'         exportPath = tempdir(),
#'         paired     = FALSE
#'     )
#' }
#'
#' @export
calculateCGEnrichment <- function(
    file = NULL, BSgenome = NULL, exportPath = NULL,
    extend = 0, shift = 0, uniq = 0,
    chr.select = NULL, paired = TRUE) {

    dataset <- eval(parse(text = paste0(BSgenome, "::", BSgenome)))

    ## Read region file
    fileName <- basename(file)
    path <- dirname(file)
    if (path == "") {
        path <- getwd()
    }
    if (!fileName %in% dir(path)) {
        stop(
            sprintf(
                "File %s not found in %s",
                shQuote(fileName), shQuote(path)
            ),
            call. = FALSE
        )
    }

    chr.lengths <- GenomeInfoDb::seqlengths(dataset)

    # uniq is applied inside the import helpers, while reads still carry their
    # real strand, as MEDIPS did. Collapsing here instead would additionally
    # merge reads sharing coordinates on opposite strands.
    if (!paired) {
        GRange.Reads <- readSingleEndFragments(
            file = file, chr.select = chr.select,
            chr.lengths = chr.lengths,
            extend = extend, shift = shift, uniq = uniq
        )
    } else {
        GRange.Reads <- readPairedFragments(
            file = file, chr.select = chr.select,
            chr.lengths = chr.lengths, uniq = uniq
        )
    }

    ## Sort chromosomes
    if (length(unique(GenomeInfoDb::seqlevels(GRange.Reads))) > 1) {
        chromosomes <- gtools::mixedsort(
            unique(GenomeInfoDb::seqlevels(GRange.Reads))
        )
    }
    if (length(unique(GenomeInfoDb::seqlevels(GRange.Reads))) == 1) {
        chromosomes <- unique(GenomeInfoDb::seqlevels(GRange.Reads))
    }

    chr_lengths <- as.numeric(GenomeInfoDb::seqlengths(dataset)[chromosomes])

    IRanges::ranges(GRange.Reads) <- IRanges::restrict(
        IRanges::ranges(GRange.Reads), +1
    )

    ## Calculate CpG density for regions
    total <- length(chromosomes)

    readsChars <- unlist(
        Biostrings::getSeq(dataset, GRange.Reads, as.character = TRUE)
    )

    # Faster to use stringr, as we are looking for exact matches.
    regions.CG <- sum(stringr::str_count(readsChars, stringr::fixed("CG")))
    regions.C <- sum(stringr::str_count(readsChars, stringr::fixed("C")))
    regions.G <- sum(stringr::str_count(readsChars, stringr::fixed("G")))
    all.genomic <- sum(stringr::str_length(readsChars))

    nReads <- length(readsChars)

    regions.relH <- as.numeric(regions.CG) / as.numeric(all.genomic) * 100
    regions.GoGe <- (
        as.numeric(regions.CG) * as.numeric(all.genomic)
    ) / (as.numeric(regions.C) * as.numeric(regions.G))

    if (BSgenome == "BSgenome.Hsapiens.NCBI.GRCh38") {
        utils::data(
            "BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution",
            package = "mesa", envir = environment()
        )
        genomicDistribution <- BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution
    } else if (BSgenome == "BSgenome.Hsapiens.UCSC.hg19") {
        utils::data(
            "BSgenome.Hsapiens.UCSC.hg19.CpG.distribution",
            package = "mesa", envir = environment()
        )
        genomicDistribution <- BSgenome.Hsapiens.UCSC.hg19.CpG.distribution
    } else {
        genomicDistribution <- calculateGenomicCGDistribution(BSgenome)
    }

    genome.relH <- genomicDistribution$genome.relH
    genome.GoGe <- genomicDistribution$genome.GoGe

    enrichment.score.relH <- regions.relH / genome.relH
    enrichment.score.GoGe <- regions.GoGe / genome.GoGe

    if (!is.null(exportPath)) {
        GRange.Reads %>%
            BiocGenerics::width() %>%
            dplyr::as_tibble() %>%
            ggplot2::ggplot(ggplot2::aes(x = value)) +
            ggplot2::geom_density(color = "black") +
            ggplot2::theme_bw() +
            ggplot2::labs(
                xlab = "Fragment Length",
                ylab = "Density",
                title = "Fragment Length Distribution",
                subtitle = stringr::str_remove(fileName, ".bam")
            )

        ggplot2::ggsave(
            file = stringr::str_replace(
                file.path(exportPath, fileName), ".bam", ".pdf"
            )
        )

        saveRDS(
            GRange.Reads,
            file = stringr::str_replace(
                file.path(exportPath, fileName), ".bam", ".rds"
            )
        )
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
    return(data.frame(
        file = file,
        relH = enrichment.score.relH,
        GoGe = enrichment.score.GoGe,
        nReads = nReads,
        nReadsWithoutPattern = numWithoutPattern,
        n100bpReads = numReads100bp,
        n100bpReadsWithoutPattern = numWithoutPatternOver100bp
    ))
}


#' Genomic positions of a motif (CG)
#'
#' Return the genomic positions of the \code{"CG"} dinucleotide for the
#' specified BSgenome and chromosomes, located directly with \pkg{Biostrings}.
#'
#' @param BSgenome Character(1). BSgenome package name.
#' @param chr.select Character vector of chromosome names to include (e.g.,
#' \code{paste0("chr", 1:22)}). If \code{NULL}, all chromosomes are used.
#'
#' @return A \link[GenomicRanges]{GRanges-class} of motif positions.
#'
#' @seealso \code{\link{calculateCGEnrichment}},
#' \code{\link{calculateCGEnrichmentGRanges}}, \pkg{Biostrings}
#'
#' @examples
#' # Requires a BSgenome package
#' # if (requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
#' #   getCGPositions("BSgenome.Hsapiens.NCBI.GRCh38", chr.select = "22")
#' # }
getCGPositions <- function(BSgenome, chr.select) {
    dataset <- eval(parse(text = paste0(BSgenome, "::", BSgenome)))

    chrs <- if (is.null(chr.select)) {
        GenomeInfoDb::seqnames(dataset)
    } else {
        as.character(chr.select)
    }

    perChr <- lapply(chrs, function(chr) {
        hits <- Biostrings::matchPattern("CG", dataset[[chr]])
        # Width 1, not 2, matching MEDIPS::MEDIPS.getPositions(), which
        # returned IRanges(start = start, end = start). Callers classify reads
        # with filter_by_non_overlaps(), where any overlap counts: a width-2
        # range would mark a read starting on the G of a CpG as containing
        # "CG" even though its extracted sequence does not.
        GenomicRanges::GRanges(
            chr,
            IRanges::IRanges(start = BiocGenerics::start(hits), width = 1L)
        )
    })

    unlist(GenomicRanges::GRangesList(perChr), use.names = FALSE)
}


#' Build a ScanBamParam, using the BAM index only when one exists
#'
#' \code{Rsamtools::scanBam()} requires an index whenever \code{which} is
#' supplied. \code{MEDIPS::getGRange()} tested for the index first and, when it
#' was absent, scanned the whole file and filtered chromosomes afterwards, so
#' unindexed BAMs stayed usable. This reproduces that behaviour.
#'
#' \code{simpleCigar = FALSE} is deliberate. \pkg{MEDIPS} defaults it to
#' \code{TRUE}, which drops reads whose CIGAR contains \code{N}, \code{S},
#' \code{H} or \code{P}, but mesa has overridden that to \code{FALSE} on both
#' \code{getGRange()} and \code{getPairedGRange()} since its first commit
#' (\code{f28d678}) - soft-clipped and spliced alignments are counted. Do not
#' "restore" the \pkg{MEDIPS} default here: it would silently drop reads.
#'
#' @param file Character(1). Path to the BAM file.
#' @param what Character vector of BAM fields to read.
#' @param flag A \code{Rsamtools::scanBamFlag()} object.
#' @param chr.select Character vector of chromosomes, or \code{NULL} for all.
#' @param chr.lengths Named numeric vector of chromosome lengths for the whole
#' genome, used as the upper bound of each scan range.
#'
#' @return A \code{list} of \code{param} (the \code{ScanBamParam}) and
#' \code{prefiltered}: \code{TRUE} when the scan is already restricted to
#' \code{chr.select}, \code{FALSE} when the caller must filter afterwards.
#'
#' @keywords internal
#' @noRd
bamScanParam <- function(file, what, flag, chr.select = NULL,
    chr.lengths = NULL) {

    hasIndex <- any(file.exists(c(
        paste0(file, ".bai"),
        sub("\\.bam$", ".bai", file, ignore.case = TRUE),
        paste0(file, ".csi")
    )))

    if (is.null(chr.select) || !hasIndex) {
        if (!is.null(chr.select)) {
            message(
                "No BAM index found for ", basename(file),
                "; scanning the whole file and selecting ",
                paste(chr.select, collapse = ", "), " afterwards."
            )
        }
        return(list(
            param = Rsamtools::ScanBamParam(
                what = what, flag = flag, simpleCigar = FALSE
            ),
            prefiltered = is.null(chr.select)
        ))
    }

    lengths <- chr.lengths[as.character(chr.select)]

    if (anyNA(lengths)) {
        stop(
            "chr.select entries absent from the BSgenome: ",
            paste(chr.select[is.na(lengths)], collapse = ", "),
            call. = FALSE
        )
    }

    which <- GenomicRanges::GRanges(
        as.character(chr.select),
        IRanges::IRanges(start = 1, end = as.integer(lengths))
    )

    list(
        param = Rsamtools::ScanBamParam(
            what = what, flag = flag, simpleCigar = FALSE, which = which
        ),
        prefiltered = TRUE
    )
}


#' Total length of the sequences being scanned
#'
#' The denominator of the \code{uniq} Poisson rate. MEDIPS used
#' \code{sum(seqlengths(dataset)[chr.select])}, which is \code{0} when
#' \code{chr.select} is \code{NULL}; fall back to the whole genome instead.
#'
#' @param chr.select Character vector of chromosomes, or \code{NULL} for all.
#' @param chr.lengths Named numeric vector of chromosome lengths for the whole
#' genome.
#'
#' @return Numeric(1).
#'
#' @keywords internal
#' @noRd
genomeLengthOf <- function(chr.select, chr.lengths) {
    if (is.null(chr.select)) {
        return(sum(as.numeric(chr.lengths)))
    }
    sum(as.numeric(chr.lengths[as.character(chr.select)]))
}


#' Apply the MEDIPS \code{uniq} duplicate-handling rule
#'
#' Reproduces the four branches of \code{MEDIPS::getGRange()} and
#' \code{MEDIPS::getPairedGRange()}.
#'
#' Must be called while \code{reads} still carries its real strand: MEDIPS
#' deduplicated before setting the strand to \code{"*"}, so collapsing
#' afterwards would additionally merge reads that share coordinates on opposite
#' strands.
#'
#' @param reads A \link[GenomicRanges]{GRanges-class} of reads.
#' @param uniq Numeric(1). \code{0} keeps every read; \code{1} keeps at most one
#' read per genomic location; a value in \code{(0, 1)} is a p-value capping the
#' reads per location at a Poisson quantile of the per-base read rate. Any
#' other value, and any logical, is an error.
#' @param genomeLength Numeric(1). Total length of the sequences scanned, the
#' denominator of the Poisson rate.
#'
#' @return A \link[GenomicRanges]{GRanges-class} of reads.
#'
#' @keywords internal
#' @noRd
dedupeReads <- function(reads, uniq, genomeLength) {

    if (is.logical(uniq)) {
        stop(
            "Parameter 'uniq' must be numeric, not logical: supply 0 to ",
            "keep all reads, 1 to keep one read per genomic location, or ",
            "a p-value in (0, 1) to cap duplicates.",
            call. = FALSE
        )
    }

    if (length(uniq) != 1 || is.na(uniq) || uniq < 0 || uniq > 1) {
        stop(
            "Parameter 'uniq' must be a single value in [0, 1]; got ",
            paste(format(uniq), collapse = ", "), ".",
            call. = FALSE
        )
    }

    if (uniq == 0) {
        return(reads)
    }

    if (uniq == 1) {
        return(BiocGenerics::unique(reads))
    }

    maxDup <- max(1, stats::qpois(1 - uniq, length(reads) / genomeLength))

    uniqReads <- BiocGenerics::unique(reads)
    dupNumber <- tabulate(
        BiocGenerics::match(reads, uniqReads),
        nbins = length(uniqReads)
    )

    rep(uniqReads, times = pmin(dupNumber, maxDup))
}


#' Import paired-end fragments from a BAM file
#'
#' Read the properly paired fragments from a BAM file and return them as a
#' \link[GenomicRanges]{GRanges-class}, one range per fragment. Replaces the
#' previous reliance on \code{MEDIPS::getPairedGRange()}, which is only usable
#' when \pkg{GenomicRanges} happens to be attached to the search path.
#'
#' Mirrors \code{MEDIPS::getPairedGRange()}: the first mate of each properly
#' mapped, non-secondary pair is taken, and the fragment span is reconstructed
#' from the leftmost of the read and its mate position plus the template length
#' (\code{isize}). \code{shift} / \code{extend} are intentionally ignored for
#' paired data (the true fragment span is known).
#'
#' @param file Character(1). Path to the BAM file. An index is used when
#' present; without one the whole file is scanned and \code{chr.select} applied
#' afterwards.
#' @param chr.select Character vector of chromosomes to import, or \code{NULL}
#' for all chromosomes.
#' @param chr.lengths Named numeric vector of chromosome lengths for the whole
#' genome, used as the upper bound of the scan range and as the denominator of
#' the \code{uniq} Poisson rate.
#' @param uniq Numeric(1). Duplicate handling, see \code{dedupeReads()}.
#'
#' @return A \link[GenomicRanges]{GRanges-class} of fragment ranges.
#'
#' @keywords internal
#' @noRd
readPairedFragments <- function(file, chr.select = NULL, chr.lengths = NULL,
    uniq = 0) {

    flag <- Rsamtools::scanBamFlag(
        isPaired = TRUE, isProperPair = TRUE,
        hasUnmappedMate = FALSE, isUnmappedQuery = FALSE,
        isFirstMateRead = TRUE, isSecondMateRead = FALSE,
        isSecondaryAlignment = FALSE
    )
    what <- c("rname", "pos", "strand", "isize", "mpos")

    scan <- bamScanParam(file, what, flag, chr.select, chr.lengths)

    readDF <- Rsamtools::scanBam(file = file, param = scan$param) %>%
        purrr::map_df(as.data.frame)

    if (!scan$prefiltered) {
        readDF <- readDF %>%
            dplyr::filter(as.character(rname) %in% as.character(chr.select))
    }

    fragments <- readDF %>%
        dplyr::mutate(
            seqnames = as.character(rname),
            start = pmin(pos, mpos),
            end = pmin(pos, mpos) + abs(isize) - 1,
            strand = as.character(strand)
        ) %>%
        plyranges::as_granges()

    fragments <- dedupeReads(
        fragments, uniq, genomeLengthOf(chr.select, chr.lengths)
    )

    BiocGenerics::strand(fragments) <- "*"

    fragments
}


#' Import single-end reads from a BAM file
#'
#' Read mapped single-end reads from a BAM file and return them as a
#' \link[GenomicRanges]{GRanges-class}, one range per read. Replaces the
#' previous reliance on \code{MEDIPS::getGRange()}, which is only usable when
#' \pkg{GenomicRanges} happens to be attached to the search path.
#'
#' Read spans are \code{[pos, pos + qwidth - 1]} (mirroring
#' \code{MEDIPS::getGRange()}). When \code{extend > 0}, reads shorter than
#' \code{extend} are lengthened to it in the 5'->3' (strand-aware) direction;
#' reads already longer are left alone.
#'
#' @param file Character(1). Path to the BAM file. An index is used when
#' present; without one the whole file is scanned and \code{chr.select} applied
#' afterwards.
#' @param chr.select Character vector of chromosomes to import, or \code{NULL}
#' for all chromosomes.
#' @param chr.lengths Named numeric vector of chromosome lengths for the whole
#' genome, used as the upper bound of the scan range and as the denominator of
#' the \code{uniq} Poisson rate.
#' @param extend Integer(1). If non-zero, reads shorter than this are extended
#' to it.
#' @param shift Integer(1). Optional strand-aware offset applied to reads.
#' @param uniq Numeric(1). Duplicate handling, see \code{dedupeReads()}.
#'
#' @return A \link[GenomicRanges]{GRanges-class} of read ranges.
#'
#' @keywords internal
#' @noRd
readSingleEndFragments <- function(file, chr.select = NULL,
    chr.lengths = NULL, extend = 0, shift = 0, uniq = 0) {

    flag <- Rsamtools::scanBamFlag(
        isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE
    )
    what <- c("rname", "pos", "strand", "qwidth")

    scan <- bamScanParam(file, what, flag, chr.select, chr.lengths)

    readDF <- Rsamtools::scanBam(file = file, param = scan$param) %>%
        purrr::map_df(as.data.frame)

    if (!scan$prefiltered) {
        readDF <- readDF %>%
            dplyr::filter(as.character(rname) %in% as.character(chr.select))
    }

    reads <- readDF %>%
        dplyr::mutate(
            seqnames = as.character(rname),
            start = pos,
            end = pos + qwidth - 1
        ) %>%
        plyranges::as_granges()

    if (shift != 0) {
        offsets <- ifelse(
            BiocGenerics::strand(reads) == "-", -shift, shift
        )
        reads <- GenomicRanges::shift(reads, offsets)
    }

    if (extend > 0) {
        # resize() on a GRanges is strand-aware: fix = "start" extends from
        # the 5' end regardless of strand (not the lower genomic
        # coordinate), so this already matches the 5'->3' extension
        # described above. pmax() reproduces the pmax(0, extend - width)
        # clamp in MEDIPS::adjustReads(): extend only ever lengthens a read,
        # it never truncates one that is already longer.
        reads <- GenomicRanges::resize(
            reads,
            width = pmax(BiocGenerics::width(reads), extend),
            fix = "start"
        )
    }

    reads <- dedupeReads(
        reads, uniq, genomeLengthOf(chr.select, chr.lengths)
    )

    BiocGenerics::strand(reads) <- "*"

    reads
}


#' CpG enrichment from GRanges of reads (MEDIPS-style)
#'
#' Compute CpG enrichment metrics (**relH** and **GoGe**) from a
#' [GenomicRanges::GRanges-class] of read spans, using the specified BSgenome.
#' Useful when reads are already represented as genomic ranges rather than
#' extracted directly from a BAM file.
#'
#' @param readGRanges `GRanges`
#' Genomic ranges representing fragments (each range = one fragment).
#'   **Default:** `NULL` (must be supplied).
#'
#' @param BSgenome `character(1)`
#' Name of a BSgenome package, e.g. `"BSgenome.Hsapiens.NCBI.GRCh38"`. The
#' package must be installed and loadable.
#'   **Default:** `NULL`.
#'
#' @param chr.select `character()` or `NULL`
#' Vector of chromosomes to restrict motif calculation.
#'   **Default:** `NULL` (use all chromosomes).
#'
#' @return A `tibble` with one row and the following columns:
#' * **relH** — sample relH normalised by genome relH.
#' * **GoGe** — sample GoGe normalised by genome GoGe.
#' * **nReads** — total reads (fragments).
#' * **nReadsWithoutPattern** — reads lacking `"CG"` motif.
#' * **n100bpReads** — reads with length >= 100 bp.
#' * **n100bpReadsWithoutPattern** — reads >= 100 bp lacking `"CG"`.
#'
#' @details
#' relH and GoGe are computed by counting `"CG"` dinucleotides in the provided
#' fragments and normalising by genome-wide expectations (from
#' [calculateGenomicCGDistribution()]). Unlike [calculateCGEnrichment()], this
#' function assumes reads are already summarised as genomic ranges.
#'
#' @seealso
#' [calculateCGEnrichment()], [calculateGenomicCGDistribution()],
#' \pkg{GenomicRanges}, \pkg{BSgenome}
#'
#' @examples
#' # Runnable toy example with synthetic reads over chr1 (requires a BSgenome)
#' if (requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)) {
#'     n <- 200
#'     gr <- GenomicRanges::GRanges(
#'         seqnames = rep(1, n),
#'         ranges = IRanges::IRanges(
#'             start = sample(1e6:2e6, n),
#'             width = sample(80:180, n, replace = TRUE)
#'         ),
#'         strand = "*"
#'     )
#'     calculateCGEnrichmentGRanges(
#'         readGRanges = gr,
#'         BSgenome    = "BSgenome.Hsapiens.NCBI.GRCh38",
#'         chr.select  = 1
#'     )
#' }
#' @export
calculateCGEnrichmentGRanges <- function(
    readGRanges = NULL, BSgenome = NULL, chr.select = NULL
) {
    dataset <- eval(parse(text = paste0(BSgenome, "::", BSgenome)))

    chromosomes <- gtools::mixedsort(
        unique(GenomeInfoDb::seqlevels(readGRanges))
    )

    chr_lengths <- as.numeric(GenomeInfoDb::seqlengths(dataset)[chromosomes])

    if (all(is.na(GenomeInfoDb::seqlengths(readGRanges)))) {
        GenomeInfoDb::seqinfo(readGRanges, pruning.mode = "coarse") <-
            GenomeInfoDb::seqinfo(BSgenome::getBSgenome(BSgenome))[
                GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(readGRanges))
            ]
        readGRanges <- IRanges::trim(readGRanges)
    }

    IRanges::ranges(readGRanges) <- IRanges::restrict(
        IRanges::ranges(readGRanges), +1
    )

    ## Calculate CpG density for regions
    total <- length(chromosomes)

    readsChars <- unlist(
        Biostrings::getSeq(dataset, readGRanges, as.character = TRUE)
    )

    # Faster to use stringr, as we are looking for exact matches.
    regions.CG <- sum(stringr::str_count(readsChars, stringr::fixed("CG")))
    regions.C <- sum(stringr::str_count(readsChars, stringr::fixed("C")))
    regions.G <- sum(stringr::str_count(readsChars, stringr::fixed("G")))
    all.genomic <- sum(stringr::str_length(readsChars))

    regions.relH <- as.numeric(regions.CG) / as.numeric(all.genomic) * 100
    regions.GoGe <- (
        as.numeric(regions.CG) * as.numeric(all.genomic)
    ) / (as.numeric(regions.C) * as.numeric(regions.G))

    if (BSgenome == "BSgenome.Hsapiens.NCBI.GRCh38") {
        utils::data(
            "BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution",
            package = "mesa", envir = environment()
        )
        genomicDistribution <- BSgenome.Hsapiens.NCBI.GRCh38.CpG.distribution
    } else if (BSgenome == "BSgenome.Mmusculus.UCSC.mm10") {
        utils::data(
            "BSgenome.Mmusculus.UCSC.mm10.CpG.distribution",
            package = "mesa", envir = environment()
        )
        genomicDistribution <- BSgenome.Mmusculus.UCSC.mm10.CpG.distribution
    } else if (BSgenome == "BSgenome.Hsapiens.UCSC.hg19") {
        utils::data(
            "BSgenome.Hsapiens.UCSC.hg19.CpG.distribution",
            package = "mesa", envir = environment()
        )
        genomicDistribution <- BSgenome.Hsapiens.UCSC.hg19.CpG.distribution
    } else {
        genomicDistribution <- calculateGenomicCGDistribution(BSgenome)
    }

    genome.relH <- genomicDistribution$genome.relH
    genome.GoGe <- genomicDistribution$genome.GoGe

    enrichment.score.relH <- regions.relH / genome.relH
    enrichment.score.GoGe <- regions.GoGe / genome.GoGe

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

    return(tibble::tibble(
        relH = enrichment.score.relH,
        GoGe = enrichment.score.GoGe,
        nReads = numReads,
        nReadsWithoutPattern = numWithoutPattern,
        n100bpReads = numReads100bp,
        n100bpReadsWithoutPattern = numWithoutPatternOver100bp
    ))
}


#' Add MEDIPS-style enrichment metrics to a qseaSet sample table
#'
#' For each sample, compute MEDIPS-style CpG enrichment metrics (**relH**,
#' **GoGe**, and read counts) from its BAM and append the results to the
#' `qseaSet` sample metadata. Supports both pulldown (default) and input
#' libraries.
#'
#' @param qseaSet `qseaSet`
#' Input object whose sample table contains BAM paths.
#'
#' @param exportPath `character(1)` or `NULL`
#' Directory to export per-sample fragment-length PDFs and read-`GRanges` RDS
#' files. If `NULL`, no files are written.
#'   **Default:** `NULL`.
#'
#' @param nonEnrich `logical(1)`
#' If `TRUE`, treat samples as **Input** libraries (columns appended under
#' `@libraries$input_file`). If `FALSE`, treat as **Pulldown** (columns appended
#' under `@libraries$file_name`).
#'   **Default:** `FALSE`.
#'
#' @param extend `integer(1)`
#' Passed to [calculateCGEnrichment()]. Extension length for unpaired reads;
#' reads already longer than this are left unchanged.
#'   **Default:** `0`.
#'
#' @param shift `integer(1)`
#' Passed to [calculateCGEnrichment()]. Shift applied to unpaired read
#' positions; unused for paired reads.
#'   **Default:** `0`.
#'
#' @param uniq `numeric(1)`
#' Passed to [calculateCGEnrichment()]. Duplicate handling: `0` keeps all
#' reads, `1` keeps one per location, a p-value in `(0, 1)` caps duplicates.
#'   **Default:** `0`.
#'
#' @param chr.select `character()` or `NULL`
#' Passed to MEDIPS range extraction; subset of chromosomes to analyse.
#'   **Default:** `NULL` (all chromosomes).
#'
#' @param paired `logical(1)`
#' Whether BAMs are paired-end (uses [MEDIPS::getPairedGRange()]).
#'   **Default:** `TRUE`.
#'
#' @param file_name `character(1)`
#' Column name in the sample table holding BAM paths when `nonEnrich = FALSE`.
#' When `nonEnrich = TRUE`, the column `input_file` is used instead.
#'   **Default:** `"file_name"`.
#'
#' @param nCores `integer(1)`
#' Number of parallel cores for `parallel::mclapply()`. Set to `1` for serial.
#'   **Default:** `1`.
#'
#' @return A `qseaSet` with new MEDIPS-style metrics appended:
#'
#' * **Sample/library metrics** (added under the appropriate library frame):
#' `relH`, `GoGe`, `nReads`, `nReadsWithoutPattern`, `n100bpReads`,
#' `n100bpReadsWithoutPattern`.
#' * **Provenance**: any export artefacts (PDF/RDS) written to `exportPath` if
#' provided.
#'
#' @details
#' Internally calls [calculateCGEnrichment()] per sample to derive relH/GoGe and
#' counts, then merges the results into the `@libraries` slot (`$file_name` for
#' pulldown, `$input_file` for input when `nonEnrich = TRUE`). Parallelisation
#' uses `parallel::mclapply()`; on non-Unix systems, `nCores` > 1 is ignored.
#'
#' @seealso
#' [calculateCGEnrichment()], [calculateCGEnrichmentGRanges()]
#'
#' @examples
#' # Requires BAM files from MEDIPSData; see \dontrun{} for a full
#' # usage example.
#' \dontrun{
#' if (requireNamespace("MEDIPSData", quietly = TRUE) &&
#'     requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
#'     bam <- system.file(
#'         "extdata",
#'         "hESCs.Input.chr22.bam",
#'         package = "MEDIPSData"
#'     )
#'
#'     data.frame(
#'         sample_name = "hESC_Input_chr22",
#'         file_name = bam,
#'         group = "Input",
#'         input_file = bam,
#'         stringsAsFactors = FALSE
#'     ) %>%
#'         qsea::createQseaSet("BSgenome.Hsapiens.UCSC.hg19") %>%
#'         addMedipsEnrichmentFactors(
#'             exportPath = tempdir(),
#'             chr.select = paste0("chr22"),
#'             paired     = FALSE
#'             # If your qsea returns a BSgenome object and as.character()
#'             # doesn't help, uncomment:
#'             # , BSgenome_pkg = "BSgenome.Hsapiens.UCSC.hg19"
#'         ) %>%
#'         qsea::getSampleTable() %>%
#'         print()
#' }
#' }
#' @export
addMedipsEnrichmentFactors <- function(
    qseaSet, exportPath = NULL, nonEnrich = FALSE,
    extend = 0, shift = 0, uniq = 0,
    chr.select = NULL, paired = TRUE,
    file_name = "file_name",
    nCores = 1) {
    BSgenome <- qsea::getParameters(qseaSet)[["BSgenome"]]

    if (nonEnrich) {
        typeString <- "non-enriched"
    } else {
        typeString <- "enriched"
    }

    message(glue::glue(
        "Adding Medips Enrichment factors to ",
        "{length(getSampleNames(qseaSet))} {typeString} samples, ",
        "using {nCores} cores."
    ))

    if (!nonEnrich) {
        colsToCheck <- c(
            "relH", "GoGe", "nReads",
            "nReadsWithoutPattern", "n100bpReadsWithoutPattern"
        )

        if (any(colsToCheck %in% colnames(qsea::getSampleTable(qseaSet)))) {
            stop(glue::glue(
                "Column {colsToCheck[",
                "colsToCheck %in% colnames(qsea::getSampleTable(qseaSet))",
                "]} already in sampleTable!"
            ))
        }
    } else {
        colsToCheck <- c(
            "input_relH", "input_GoGe", "input_nReads",
            "input_nReadsWithoutPattern", "input_n100bpReadsWithoutPattern"
        )

        if (any(colsToCheck %in% colnames(qsea::getSampleTable(qseaSet)))) {
            stop(glue::glue(
                "Column {colsToCheck[",
                "colsToCheck %in% colnames(qsea::getSampleTable(qseaSet))",
                "]} already in sampleTable!"
            ))
        }
    }

    if (!nonEnrich) {
        fileNames <- qsea::getSampleTable(qseaSet) %>% dplyr::pull(file_name)
    } else {
        fileNames <- qsea::getSampleTable(qseaSet) %>% dplyr::pull(input_file)
    }


    enrichData <- parallel::mclapply(fileNames,
        function(x) {
            calculateCGEnrichment(x,
                BSgenome = BSgenome, exportPath = exportPath,
                extend = extend, shift = shift, uniq = uniq,
                chr.select = chr.select, paired = paired
            )
        },
        mc.cores = nCores
    ) %>%
        do.call(rbind, .)

    if (!nonEnrich) {
        qseaSet@libraries$file_name <- qseaSet@libraries$file_name %>%
            cbind(dplyr::select(enrichData, -file))
    } else {
        qseaSet@libraries$input_file <- qseaSet@libraries$input_file %>%
            cbind(dplyr::select(enrichData, -file))
    }

    return(qseaSet)
}
