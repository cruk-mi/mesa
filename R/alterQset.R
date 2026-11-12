#' Subset a qseaSet by overlaps / non-overlaps with genomic regions
#'
#' `filterByOverlaps()` retains windows that **overlap**, while
#' `filterByNonOverlaps()` retains windows that **do not overlap** a set of
#' genomic regions. Inputs may use different chromosome styles (with or without
#' the `"chr"` prefix); styles are harmonized automatically before overlap
#' computation.
#'
#' @param qseaSet `qseaSet`
#' A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param regionsToOverlap `GRanges` or `data.frame`
#' Genomic regions used for filtering. If a `data.frame`, it must contain
#' columns `seqnames` (`character()`), `start` (`integer()`), and `end`
#' (`integer()`); the object is coerced internally to a `GRanges`.
#'
#' @return A `qseaSet` object:
#'
#' * **filterByOverlaps()**: returns the input `qseaSet` restricted to
#' windows that **overlap** `regionsToOverlap`.
#'
#' * **filterByNonOverlaps()**: returns the input `qseaSet` restricted to
#' windows that **do not overlap** `regionsToOverlap`.
#'
#' @details
#' Chromosome naming is aligned using [`GenomeInfoDb::seqlevelsStyle()`] so that
#' overlap computation is consistent even when one input uses `"chr1"` and the
#' other uses `"1"`. Only the **window set** is filtered; sample metadata and
#' counts are subset accordingly.
#'
#' @seealso
#' \code{\link[plyranges]{filter_by_overlaps}},
#' \code{\link[plyranges]{filter_by_non_overlaps}},
#' \code{\link[GenomicRanges]{GRanges-class}}
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Data frame input
#' rgns <- data.frame(seqnames = 7L, start = 25002001L, end = 25017900L)
#' filterByOverlaps(exampleTumourNormal, regionsToOverlap = rgns)
#' filterByNonOverlaps(exampleTumourNormal, regionsToOverlap = rgns)
#'
#' # GRanges input
#' gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(25002001L, 25017900L))
#' filterByOverlaps(exampleTumourNormal, regionsToOverlap = gr)
#' filterByNonOverlaps(exampleTumourNormal, regionsToOverlap = gr)
#'
#' @rdname alterQsetOverlap
#' @export
filterByOverlaps <- function(qseaSet, regionsToOverlap) {
    if (is.data.frame(regionsToOverlap)) {
        if (length(intersect(
            colnames(regionsToOverlap), c("seqnames", "start", "end")
        )) != 3) {
            stop(
                "regionsToOverlap must be a GRanges object",
                " or a dataframe with seqnames, start and end."
            )
        }
    }

    regionsToOverlap <- plyranges::as_granges(regionsToOverlap)

    qseaSetChr <- qseaSet %>%
        qsea::getRegions() %>%
        GenomeInfoDb::seqinfo() %>%
        GenomeInfoDb::seqnames() %>%
        stringr::str_detect("chr") %>%
        any()
    windowsChr <- regionsToOverlap %>%
        GenomeInfoDb::seqinfo() %>%
        GenomeInfoDb::seqnames() %>%
        stringr::str_detect("chr") %>%
        any()

    if (qseaSetChr & !windowsChr) {
        regionsToOverlap <- regionsToOverlap %>%
            tibble::as_tibble() %>%
            dplyr::mutate(seqnames = paste0("chr", seqnames)) %>%
            plyranges::as_granges()
    }

    if (!qseaSetChr & windowsChr) {
        regionsToOverlap <- regionsToOverlap %>%
            tibble::as_tibble() %>%
            dplyr::mutate(seqnames = stringr::str_remove(seqnames, "chr")) %>%
            plyranges::as_granges()
    }

    windowsToKeep <- which(
        plyranges::count_overlaps(
            qsea::getRegions(qseaSet), regionsToOverlap
        ) > 0
    )

    qseaSet <- qseaSet %>%
        setCounts(
            qsea::getCounts(qseaSet)[windowsToKeep, , drop = FALSE]
        ) %>%
        setRegions(
            plyranges::filter_by_overlaps(
                qsea::getRegions(qseaSet), regionsToOverlap
            )
        )

    return(qseaSet)
}


#' @rdname alterQsetOverlap
#' @export
filterByNonOverlaps <- function(qseaSet, regionsToOverlap) {
    if (is.data.frame(regionsToOverlap)) {
        if (length(intersect(
            colnames(regionsToOverlap), c("seqnames", "start", "end")
        )) != 3) {
            stop(
                "regionsToOverlap must be a GRanges object or a",
                " dataframe with seqnames, start and end."
            )
        }
    }

    regionsToOverlap <- plyranges::as_granges(regionsToOverlap)

    qseaSetChr <- qseaSet %>%
        qsea::getRegions() %>%
        GenomeInfoDb::seqinfo() %>%
        GenomeInfoDb::seqnames() %>%
        stringr::str_detect("chr") %>%
        any()
    windowsChr <- regionsToOverlap %>%
        GenomeInfoDb::seqinfo() %>%
        GenomeInfoDb::seqnames() %>%
        stringr::str_detect("chr") %>%
        any()

    if (qseaSetChr & !windowsChr) {
        regionsToOverlap <- regionsToOverlap %>%
            tibble::as_tibble() %>%
            dplyr::mutate(seqnames = paste0("chr", seqnames)) %>%
            plyranges::as_granges()
    }

    if (!qseaSetChr & windowsChr) {
        regionsToOverlap <- regionsToOverlap %>%
            tibble::as_tibble() %>%
            dplyr::mutate(seqnames = stringr::str_remove(seqnames, "chr")) %>%
            plyranges::as_granges()
    }

    GRangesToKeep <- qsea::getRegions(qseaSet) %>%
        plyranges::filter_by_non_overlaps(
            plyranges::as_granges(regionsToOverlap)
        )

    qseaSet %>%
        filterByOverlaps(GRangesToKeep) %>%
        return()
}


#' Add library metrics to the qseaSet sample table
#'
#' Append library information to the `sampleTable` of a `qseaSet`. Columns are
#' taken from the `"file_name"` table of [getLibraries()], and when available,
#' from the `"input_file"` table (prefixed with `"input_"`). Existing columns
#' are not duplicated.
#'
#' @param qseaSet `qseaSet`
#' A qseaSet object containing sequencing data.
#'
#' @return A `qseaSet` object:
#'
#' * **addLibraryInformation()**: returns the same `qseaSet` with its
#' `sampleTable` augmented by library-level columns.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' qs <- addLibraryInformation(exampleTumourNormal)
#' head(qsea::getSampleTable(qs))
#'
#' @export
addLibraryInformation <- function(qseaSet) {
    curColNames <- qseaSet %>%
        qsea::getSampleTable() %>%
        colnames()

    qseaSet <- setSampleTable(qseaSet, cbind(
        qsea::getSampleTable(qseaSet),
        getLibraries(qseaSet, "file_name") %>%
            dplyr::select(-tidyselect::any_of(curColNames))
    ))

    if ("input_file" %in% names(getLibraries(qseaSet))) {
        curColNames <- qseaSet %>%
            qsea::getSampleTable() %>%
            colnames()

        qseaSet <- setSampleTable(qseaSet, cbind(
            qsea::getSampleTable(qseaSet),
            getLibraries(qseaSet, "input_file") %>%
                dplyr::rename_with(~ paste0("input_", .x)) %>%
                dplyr::select(-tidyselect::any_of(curColNames)) %>%
                dplyr::select(-tidyselect::matches("^library_factor$|^offset$"))
        ))
    }


    return(qseaSet)
}


#' Subset a qseaSet by samples
#'
#' Restrict a `qseaSet` to a specified set of samples. All internal matrices and
#' slots are subset consistently. Use exactly one of `samplesToKeep` or
#' `samplesToDrop`.
#'
#' @param qseaSet `qseaSet`
#' A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param samplesToKeep `character()`
#' Names of samples to retain. Cannot be used together with `samplesToDrop`.
#'   **Default:** `NULL`.
#'
#' @param samplesToDrop `character()`
#' Names of samples to remove. Cannot be used together with `samplesToKeep`.
#'   **Default:** `NULL`.
#'
#' @return A `qseaSet` object:
#'
#' * **subsetQset()**: returns the input `qseaSet` restricted to the
#' selected samples.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Keep only two samples
#' subsetQset(exampleTumourNormal, samplesToKeep = c("Colon1_T", "Colon1_N"))
#'
#' # Drop one sample
#' subsetQset(exampleTumourNormal, samplesToDrop = "Lung1_N")
#'
#' @export
subsetQset <- function(qseaSet, samplesToKeep = NULL, samplesToDrop = NULL) {
    if (length(samplesToKeep) == 0 & length(samplesToDrop) == 0) {
        message("No samples remaining, returning an empty qseaSet.")
    }

    if (length(samplesToKeep) > 0 & length(samplesToDrop) > 0) {
        stop("Can only specify samples to keep or to drop, not both.")
    }

    samplesNotInSet <- setdiff(
        c(samplesToKeep, samplesToDrop), qsea::getSampleNames(qseaSet)
    )
    if (length(samplesNotInSet) > 0) {
        stop(glue::glue("Sample {samplesNotInSet} not present in the qseaSet!

                    "))
    }

    if (length(samplesToDrop) > 0) {
        samplesToKeep <- setdiff(qsea::getSampleNames(qseaSet), samplesToDrop)
    }

    newSet <- qseaSet %>%
        setSampleTable(
            qsea::getSampleTable(qseaSet)[samplesToKeep, , drop = FALSE]
        ) %>%
        setCounts(qsea::getCounts(qseaSet)[, samplesToKeep, drop = FALSE]) %>%
        qsea::setZygosity(
            qsea::getZygosity(qseaSet)[samplesToKeep, , drop = FALSE]
        ) %>%
        setLibraries(
            getLibraries(qseaSet, "file_name")[samplesToKeep, , drop = FALSE],
            "file_name"
        ) %>%
        setLibraries(
            getLibraries(qseaSet, "input_file")[samplesToKeep, , drop = FALSE],
            "input_file"
        ) %>%
        setEnrichment(
            getEnrichment(qseaSet, "parameters")[
                samplesToKeep, ,
                drop = FALSE
            ],
            "parameters"
        ) %>%
        setEnrichment(
            getEnrichment(qseaSet, "factors")[, samplesToKeep, drop = FALSE],
            "factors"
        )

    if (qsea::hasCNV(qseaSet)) {
        newSet <- setCNV(newSet, qsea::getCNV(qseaSet, samplesToKeep))
    }

    return(newSet)
}


#' Rename samples in a qseaSet by regex
#'
#' Apply a regular expression replacement to sample names and update all
#' relevant slots consistently (counts, CNV, enrichment, zygosity, libraries,
#' and the `sampleTable`).
#'
#' @param qseaSet `qseaSet`
#' A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param pattern `character(1)`
#' Regular expression used to match substrings in sample names.
#'
#' @param replacement `character(1)`
#' String to replace matches of `pattern`. Defaults to `""` (empty string).
#'
#' @return A `qseaSet` object:
#'
#' * **renameQsetNames()**: returns the input `qseaSet` with sample names
#' updated according to `pattern` and `replacement`. All dependent slots are
#' updated consistently.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Replace "T" with "Tumour" in sample names
#' renameQsetNames(exampleTumourNormal,
#'     pattern = "T",
#'     replacement = "Tumour"
#' )
#'
#' @export
renameQsetNames <- function(qseaSet, pattern, replacement = "") {
    renamedNames <- qseaSet %>%
        qsea::getSampleTable() %>%
        dplyr::pull(sample_name) %>%
        stringr::str_replace_all(pattern, replacement)

    if (all(renamedNames == qsea::getSampleNames(qseaSet))) {
        message("Renaming had no effect!")
        return(qseaSet)
    }

    asValidNames <- base::make.names(renamedNames)

    if (any(asValidNames != renamedNames)) {
        invalidNames <- glue::glue_collapse(
            renamedNames[renamedNames != asValidNames],
            sep = "\n    "
        )
        stop(glue::glue(
            "Sample names must be valid names for columns in R without ",
            "quoting.\n",
            "See the help for base::make.names, but generally use only ",
            "letters, numbers,\n",
            "underscores and dots, and names can't start with a number.\n",
            "Issues were found with:\n    {invalidNames}"
        ))
    }

    if (any(duplicated(renamedNames))) {
        stop(glue::glue(
            "Duplicate sample_name now present: ",
            "{renamedNames[duplicated(renamedNames)]}"
        ))
    }

    return(.setSampleNames(qseaSet, renamedNames))
}


#' Pool (merge) samples that share a common name prefix
#'
#' Merge samples whose names share a common prefix by removing a suffix pattern
#' (`mergeString`) and summing counts across the resulting groups. CNV and
#' zygosity values are combined using fragment‐count weights. The result is a
#' new `qseaSet` with pooled samples and updated metadata.
#'
#' @param qseaSet `qseaSet`
#' A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param mergeString `character(1)`
#' Regular expression used to remove the **suffix** that distinguishes samples
#' to be merged. The remaining prefix becomes the pooled sample name, and the
#' samples that give the same prefix are pooled together. For example, with
#' names like `"Patient1_T"` and `"Patient1_N"`, use `mergeString = "_[TN]$"` to
#' pool them into `"Patient1"`.
#'
#' @param renormalise `logical(1)`
#' Whether to run [addNormalisation()] on the result. The pooled samples have no
#' library factors, offsets or enrichment parameters until it has run, so beta
#' values and DMR calls need it. Set this to `FALSE` only for a set that cannot
#' be normalised, such as one whose windows are all CpG-dense.
#'   **Default:** `TRUE`.
#'
#' @param maxPatternDensity `numeric(1)`
#' Maximum pattern density in a window to consider it for the background
#' calculation. Passed to [addNormalisation()] when `renormalise = TRUE`.
#'   **Default:** `0.05`.
#'
#' @return A `qseaSet` object:
#'
#' * **poolSamples()**: returns a new `qseaSet` where samples with names sharing
#' a common prefix (after applying `mergeString`) are merged. Counts are summed,
#' CNV and zygosity values are combined using weighted averages, and library
#' statistics are aggregated.
#'
#' @details
#' Internally, pooled counts are simple sums. CNV values and zygosity are
#' combined by a weighted average using `total_fragments` from the
#' `"input_file"` table of [getLibraries()], or from the `"file_name"` table
#' when the set has no input libraries. Library statistics are aggregated
#' (sums for counts, weighted means for rates). When `renormalise = TRUE`,
#' [addNormalisation()] is run with the `enrichmentMethod` recorded in the
#' input set (`"blind1-15"` if there is none).
#'
#' @seealso [addNormalisation()], [mixSamples()]
#'
#' @examples
#' data(exampleTumourNormal)
#' # Pool Tumour/Normal pairs by patient into a single sample per patient:
#' # e.g., "Colon1_T","Colon1_N" -> "Colon1". The windows of the example data
#' # are too CpG-dense to estimate the background, so skip the normalisation.
#' poolSamples(exampleTumourNormal, mergeString = "_[TN]$", renormalise = FALSE)
#'
#' # A set with a range of CpG densities can be renormalised after pooling
#' qsea::getExampleQseaSet(expSamplingDepth = 1e5) %>%
#'     poolSamples(mergeString = "[TN]$", maxPatternDensity = 0.5)
#' @export
poolSamples <- function(qseaSet, mergeString, renormalise = TRUE,
    maxPatternDensity = 0.05) {
    ## TODO: rewrite this function to use a column.
    ## TODO: make sure it works in general

    samplesToMerge <- qsea::getSampleTable(qseaSet) %>%
        dplyr::pull(sample_name) %>%
        stringr::str_subset(mergeString)

    poolNames <- stringr::str_remove(samplesToMerge, mergeString)
    newNames <- unique(poolNames)

    if (length(newNames) >= length(samplesToMerge)) {
        stop(glue::glue("New names are the same length as the old ones!"))
    }

    # the samples that make up each pool, matched exactly
    poolMembers <- split(samplesToMerge, factor(poolNames, levels = newNames))

    newSet <- qseaSet %>%
        subsetQset(samplesToKeep = samplesToMerge)

    counts <- qsea::getCounts(qseaSet)
    zygosity <- qsea::getZygosity(qseaSet)

    newCounts <- vapply(
        newNames,
        function(x) {
            rowSums(counts[, poolMembers[[x]], drop = FALSE])
        },
        FUN.VALUE = numeric(nrow(counts))
    )

    # weight CNV and zygosity by the input libraries, or by the pulldown
    # libraries when there are no inputs
    weightLibraries <- getLibraries(qseaSet, "input_file")
    if (is.null(weightLibraries)) {
        weightLibraries <- getLibraries(qseaSet, "file_name")
    }

    poolWeights <- lapply(poolMembers, function(members) {
        weights <- weightLibraries[members, "total_fragments"]
        weights / sum(weights)
    })

    cnvMat <- qseaSet %>%
        qsea::getCNV() %>%
        GenomicRanges::mcols() %>%
        as.matrix()

    newCNVvals <- vapply(
        newNames,
        function(x) {
            weights <- poolWeights[[x]]
            rowSums(
                cnvMat[, poolMembers[[x]], drop = FALSE] %*%
                    diag(weights, nrow = length(weights))
            )
        },
        FUN.VALUE = numeric(nrow(cnvMat))
    )

    newZygosity <- vapply(
        newNames,
        function(x) {
            colSums(
                zygosity[poolMembers[[x]], , drop = FALSE] *
                    poolWeights[[x]]
            )
        },
        FUN.VALUE = numeric(ncol(zygosity))
    )
    # vapply() drops to a vector when there is one chromosome
    newZygosity <- t(matrix(
        newZygosity,
        ncol = length(newNames),
        dimnames = list(colnames(zygosity), newNames)
    ))

    newCNV <- qsea::getCNV(newSet)
    GenomicRanges::mcols(newCNV) <- newCNVvals

    newSet <- newSet %>%
        setCounts(newCounts) %>%
        setCNV(newCNV) %>%
        qsea::setZygosity(newZygosity)

    newSampleTable <- purrr::map_dfr(newNames, function(x) {
        qsea::getSampleTable(qseaSet)[poolMembers[[x]], , drop = FALSE] %>%
            dplyr::select(-tidyselect::matches("file_name|input_file")) %>%
            dplyr::mutate(
                sample_name = x,
                rownameCol = x
            ) %>%
            dplyr::distinct() %>%
            dplyr::slice(1) %>%
            tibble::remove_rownames() %>%
            tibble::column_to_rownames("rownameCol")
    })

    newSet <- setSampleTable(newSet, newSampleTable)

    poolLibraries <- function(libraries) {
        purrr::map_dfr(newNames, function(x) {
            libraries[poolMembers[[x]], , drop = FALSE] %>%
                tibble::as_tibble() %>%
                dplyr::summarise(
                    dplyr::across(
                        tidyselect::where(is.numeric) &
                            !tidyselect::matches("reads|pairs|r1s|fragments"),
                        ~ stats::weighted.mean(.x, total_fragments)
                    ),
                    dplyr::across(
                        tidyselect::matches("reads|pairs|r1s|fragments"),
                        ~ sum(.x)
                    )
                ) %>%
                dplyr::mutate(
                    rownameCol = x, library_factor = NA, offset = NA
                ) %>%
                tibble::remove_rownames() %>%
                tibble::column_to_rownames("rownameCol")
        }) %>%
            as.data.frame()
    }

    newSet <- setLibraries(
        newSet, poolLibraries(getLibraries(qseaSet, "file_name")), "file_name"
    )
    if (!is.null(getLibraries(qseaSet, "input_file"))) {
        newSet <- setLibraries(
            newSet,
            poolLibraries(getLibraries(qseaSet, "input_file")),
            "input_file"
        )
    }

    # the pooled samples have no enrichment fit until addNormalisation() runs
    enrichParams <- getEnrichment(qseaSet, "parameters")
    if (!is.null(enrichParams)) {
        enrichFactors <- getEnrichment(qseaSet, "factors")
        newSet <- newSet %>%
            setEnrichment(
                matrix(
                    NA_real_,
                    nrow = length(newNames), ncol = ncol(enrichParams),
                    dimnames = list(newNames, colnames(enrichParams))
                ),
                "parameters"
            ) %>%
            setEnrichment(
                matrix(
                    NA_real_,
                    nrow = nrow(enrichFactors), ncol = length(newNames),
                    dimnames = list(rownames(enrichFactors), newNames)
                ),
                "factors"
            )
    }

    oldSet <- qseaSet %>%
        subsetQset(samplesToDrop = samplesToMerge)

    finalSet <- combineQsets(oldSet, newSet)

    if (renormalise) {
        enrichmentMethod <- qsea::getParameters(qseaSet)$enrichmentMethod
        if (is.null(enrichmentMethod)) {
            enrichmentMethod <- "blind1-15"
        }
        finalSet <- addNormalisation(
            finalSet,
            enrichmentMethod = enrichmentMethod,
            maxPatternDensity = maxPatternDensity
        )
    }

    return(finalSet)
}


#' Rename samples using a column from the sample table
#'
#' Replace sample names with values from a user-supplied column in the
#' `sampleTable` of a `qseaSet`, updating all dependent slots consistently.
#'
#' @param qseaSet `qseaSet`
#' A qseaSet object containing methylation-enriched sequencing data.
#'
#' @param newNameColumn `character(1)`
#' The name of a column in the `sampleTable` containing new sample names. Values
#' must be unique and contain no `NA`s. (Unquoted column names are supported via
#' tidy evaluation.)
#'
#' @return A `qseaSet` object:
#'
#' * **renameSamples()**: returns the input `qseaSet` with sample names replaced
#' by values from `newNameColumn`. All associated slots are updated to reflect
#' the new naming.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' exampleTumourNormal %>%
#'     dplyr::mutate(
#'         newName = paste0(
#'             "Sample",
#'             seq_len(nrow(qsea::getSampleTable(.)))
#'         )
#'     ) %>%
#'     renameSamples(newNameColumn = "newName")
#'
#' @export
renameSamples <- function(qseaSet, newNameColumn) {
    newNameColumn <- rlang::enquo(newNameColumn)

    renamedNames <- qseaSet %>%
        qsea::getSampleTable() %>%
        dplyr::pull(!!newNameColumn)

    if (any(is.na(renamedNames))) {
        stop(glue::glue("NAs present in the new sample name column"))
    }

    asValidNames <- base::make.names(renamedNames)

    if (any(asValidNames != renamedNames)) {
        invalidNames <- glue::glue_collapse(
            renamedNames[renamedNames != asValidNames],
            sep = "\n    "
        )
        stop(glue::glue(
            "Sample names must be valid names for columns in R without ",
            "quoting.\n",
            "See the help for base::make.names, but generally use only ",
            "letters, numbers,\n",
            "underscores and dots, and names can't start with a number.\n",
            "Issues were found with:\n    {invalidNames}"
        ))
    }

    if (all(renamedNames == qsea::getSampleNames(qseaSet))) {
        message("Renaming had no effect!")
        return(qseaSet)
    }
    if (any(duplicated(renamedNames))) {
        stop(glue::glue(
            "Duplicate sample_name now present: ",
            "{renamedNames[duplicated(renamedNames)]}"
        ))
    }
    return(.setSampleNames(qseaSet, renamedNames))
}
