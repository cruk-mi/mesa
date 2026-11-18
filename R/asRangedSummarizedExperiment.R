#' Convert a qseaSet to a RangedSummarizedExperiment
#'
#' Creates a [SummarizedExperiment::RangedSummarizedExperiment-class] from a
#' `qseaSet`, for use with the wider Bioconductor ecosystem. The conversion is
#' one-way and not guaranteed to be lossless: the `qseaSet` remains the primary
#' object in mesa.
#'
#' The slots map as follows:
#'
#' * windows -> `rowRanges`
#' * counts, NRPM, beta -> `assays` (`counts`, `nrpm`, `beta`)
#' * CNV -> `assays$cnv`, expanded from the CNV bins to every window
#' * sample table, including the per-sample library information added by
#'   [addLibraryInformation()] (which is called if it has not been already),
#'   and the per-sample enrichment parameters (columns `enrichment.<name>`)
#'   -> `colData`
#' * parameters, the CpG-density enrichment model (`factors`, `density`, `n`)
#'   and the original CNV bins -> `metadata()`
#'
#' @param qs A `qseaSet`.
#' @param assays Which assays to include, any of `"counts"`, `"nrpm"`,
#'   `"beta"`, `"cnv"`. Defaults to `"counts"` only, as genome-wide
#'   objects can be very large.
#'
#' @return A `RangedSummarizedExperiment` with windows as rows and samples as
#'   columns.
#'
#' @examples
#' se <- asRangedSummarizedExperiment(exampleTumourNormal)
#' se
#'
#' se <- asRangedSummarizedExperiment(
#'     exampleTumourNormal,
#'     assays = c("counts", "nrpm", "beta", "cnv")
#' )
#' SummarizedExperiment::assayNames(se)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges findOverlaps pintersect
#' @importFrom IRanges width
#' @importFrom S4Vectors queryHits subjectHits
#' @export
asRangedSummarizedExperiment <- function(qs, assays = "counts") {
    if (!methods::is(qs, "qseaSet")) {
        stop("`qs` must be a qseaSet.", call. = FALSE)
    }
    assays <- match.arg(assays, c("counts", "nrpm", "beta", "cnv"),
                        several.ok = TRUE)

    samples <- getSampleNames(qs)
    windows <- getRegions(qs)
    winNames <- getWindowNames(qs)
    if (length(winNames) == length(windows)) {
        names(windows) <- winNames
    }

    assayList <- list()
    if ("counts" %in% assays) {
        assayList$counts <- .tableToMatrix(getCountTable(qs), samples)
    }
    if ("nrpm" %in% assays) {
        assayList$nrpm <- .tableToMatrix(getNRPMTable(qs), samples)
    }
    if ("beta" %in% assays) {
        assayList$beta <- .tableToMatrix(getBetaTable(qs), samples)
    }

    cnv <- getCNV(qs)
    if ("cnv" %in% assays) {
        if (is.null(cnv) || length(cnv) == 0) {
            warning("No CNV found in the qseaSet, skipping the cnv assay.",
                    call. = FALSE)
        } else {
            assayList$cnv <- .expandCNV(cnv, windows, samples)
        }
    }

    assayList <- lapply(assayList, function(m) {
        rownames(m) <- names(windows)
        m
    })

    # Per-sample library QC is moved onto the sample table by mesa's own
    # addLibraryInformation(), so the column names match the rest of mesa.
    if (!"total_fragments" %in% colnames(getSampleTable(qs))) {
        qs <- addLibraryInformation(qs)
    }

    colData <- as.data.frame(getSampleTable(qs))
    rownames(colData) <- samples

    enrichment <- .tryGet(getEnrichment, qs)

    # Per-sample enrichment parameters also go into colData, without
    # overwriting any column that is already there
    extra <- .enrichmentColumns(enrichment, samples)
    extra <- extra[!names(extra) %in% colnames(colData)]
    for (nm in names(extra)) {
        colData[[nm]] <- extra[[nm]]
    }

    # Bin-level enrichment model (indexed by CpG density, not by window or
    # sample), plus the original CNV bins, go into metadata
    meta <- list(
        parameters = getParameters(qs),
        enrichment = enrichment,
        cnvBins = cnv
    )
    meta <- meta[!vapply(meta, is.null, logical(1))]

    SummarizedExperiment::SummarizedExperiment(
        assays = assayList,
        rowRanges = windows,
        colData = colData,
        metadata = meta
    )
}

# Turn a mesa data table (a matrix, or a wide table with extra coordinate
# columns) into a numeric matrix with one column per sample, in sample order.
.tableToMatrix <- function(tab, samples) {
    if (is.matrix(tab)) {
        m <- tab
    } else {
        tab <- as.data.frame(tab)
        missingSamples <- setdiff(samples, colnames(tab))
        if (length(missingSamples) > 0) {
            stop(
                "Sample columns not found in table: ",
                paste(missingSamples, collapse = ", "), call. = FALSE
            )
        }
        m <- as.matrix(tab[, samples, drop = FALSE])
    }
    m[, samples, drop = FALSE]
}

# Copy each CNV bin's value to the windows it covers. If a window overlaps
# several bins, take the one with the largest overlap. Uncovered windows
# are NA.
.expandCNV <- function(cnv, windows, samples) {
    cnvVals <- as.matrix(as.data.frame(S4Vectors::mcols(cnv)))
    missingSamples <- setdiff(samples, colnames(cnvVals))
    if (length(missingSamples) > 0) {
        stop(
            "CNV has no values for samples: ",
            paste(missingSamples, collapse = ", "), call. = FALSE
        )
    }
    cnvVals <- cnvVals[, samples, drop = FALSE]

    hits <- GenomicRanges::findOverlaps(windows, cnv, ignore.strand = TRUE)
    q <- S4Vectors::queryHits(hits)
    s <- S4Vectors::subjectHits(hits)
    ov <- IRanges::width(
        GenomicRanges::pintersect(windows[q], cnv[s], ignore.strand = TRUE)
    )

    ord <- order(q, -ov)
    q <- q[ord]
    s <- s[ord]
    keep <- !duplicated(q)

    out <- matrix(
        NA_real_, nrow = length(windows), ncol = length(samples),
        dimnames = list(NULL, samples)
    )
    out[q[keep], ] <- cnvVals[s[keep], , drop = FALSE]
    out
}

# The per-sample enrichment parameters (one row per sample).
.enrichmentColumns <- function(enrichment, samples) {
    out <- list()
    pars <- enrichment$parameters
    if (is.null(pars) || NROW(pars) != length(samples)) {
        return(out)
    }
    pars <- as.matrix(pars)
    if (!is.null(rownames(pars)) && all(samples %in% rownames(pars))) {
        pars <- pars[samples, , drop = FALSE]
    }
    parNames <- colnames(pars)
    if (is.null(parNames)) {
        parNames <- seq_len(ncol(pars))
    }
    for (i in seq_len(ncol(pars))) {
        out[[paste("enrichment", parNames[i], sep = ".")]] <- pars[, i]
    }
    out
}

.tryGet <- function(getter, qs) {
    tryCatch(getter(qs), error = function(e) NULL)
}