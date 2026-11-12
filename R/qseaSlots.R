# The only place in mesa that reaches into a qseaSet slot with `@`. Everything
# else goes through these accessors or through the qsea::get* ones.

#' Get and set the slots of a qseaSet
#'
#' Read and replace the contents of a `qseaSet`. qsea exports accessors for
#' most slots, but none for `libraries` or `enrichment`, and no replacement
#' functions at all. These generics fill those gaps, so that code never has to
#' reach into a slot with `@`.
#'
#' @name qseaSlots
#'
#' @param object `qseaSet`
#'   The object to query or modify.
#'
#' @param type `character(1)` or `NULL`
#'   Which library table to return or replace: `"file_name"` for the pulldown
#'   libraries, `"input_file"` for the input libraries, or `"input_counts"`
#'   for the per-window input counts. `NULL` returns the whole list.
#'   **Default:** `NULL`.
#'
#' @param name `character(1)` or `NULL`
#'   Which element of the `enrichment` slot to return or replace, such as
#'   `"parameters"`, `"factors"`, `"pattern_name"`, `"density"` or `"n"`.
#'   `NULL` returns or replaces the whole list. **Default:** `NULL`.
#'
#' @param value
#'   The replacement value: a `matrix` for `setCounts()`, a `GRanges` for
#'   `setRegions()` and `setCNV()`, a `data.frame` for `setSampleTable()`, and
#'   a `data.frame`/`matrix` (or a `list` when `type`/`name` is `NULL`) for
#'   `setLibraries()` and `setEnrichment()`.
#'
#' @param values `list`
#'   Named list of parameters to merge into the `parameters` slot. Entries
#'   with an existing name are overwritten, the rest are left alone.
#'
#' @param ...
#'   Additional arguments passed to methods. Not used by the `qseaSet`
#'   methods.
#'
#' @return
#' * **getLibraries()**: the `libraries` slot, a named `list`, or the single
#'   table named by `type`. Returns `NULL` when the set has no such table.
#' * **getEnrichment()**: the `enrichment` slot, a named `list`, or the
#'   single element named by `name`. Returns `NULL` when absent.
#' * **setCounts()**, **setRegions()**, **setSampleTable()**,
#'   **setParameters()**, **setLibraries()**, **setEnrichment()** and
#'   **setCNV()**: the updated `qseaSet`.
#'
#' @details
#' The setters replace a slot without checking it against the rest of the
#' object, so a sequence of calls can leave a `qseaSet` temporarily
#' inconsistent (different numbers of samples in the count matrix and the
#' sample table, say). Replace every affected slot before using the result.
#'
#' Zygosity has no setter here because qsea already exports one,
#' [qsea::setZygosity()].
#'
#' @seealso
#' [qsea::getCounts()], [qsea::getRegions()], [qsea::getSampleTable()],
#' [qsea::getParameters()], [qsea::getCNV()], [qsea::getZygosity()],
#' [qsea::setZygosity()], [subsetQset()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # The pulldown library table, and the whole libraries slot
#' head(getLibraries(exampleTumourNormal, "file_name"))
#' names(getLibraries(exampleTumourNormal))
#'
#' # The fitted enrichment parameters
#' head(getEnrichment(exampleTumourNormal, "parameters"))
#'
#' # Replacing a slot returns a new qseaSet
#' sampleTable <- qsea::getSampleTable(exampleTumourNormal)
#' sampleTable$batch <- "A"
#' qs <- setSampleTable(exampleTumourNormal, sampleTable)
#' head(qsea::getSampleTable(qs))
#'
#' # Parameters are merged, not replaced wholesale
#' qs <- setParameters(exampleTumourNormal, list(minMapQual = 10))
#' qsea::getParameters(qs)$minMapQual
NULL


#' @rdname qseaSlots
#' @export
setGeneric("getLibraries", function(object, ...) {
    standardGeneric("getLibraries")
})

#' @rdname qseaSlots
#' @export
setMethod("getLibraries", "qseaSet", function(object, type = NULL) {
    if (is.null(type)) {
        return(object@libraries)
    }
    object@libraries[[type]]
})


#' @rdname qseaSlots
#' @export
setGeneric("setLibraries", function(object, ...) {
    standardGeneric("setLibraries")
})

#' @rdname qseaSlots
#' @export
setMethod("setLibraries", "qseaSet", function(object, value, type = NULL) {
    if (is.null(type)) {
        object@libraries <- value
    } else {
        object@libraries[[type]] <- value
    }
    object
})


#' @rdname qseaSlots
#' @export
setGeneric("getEnrichment", function(object, ...) {
    standardGeneric("getEnrichment")
})

#' @rdname qseaSlots
#' @export
setMethod("getEnrichment", "qseaSet", function(object, name = NULL) {
    if (is.null(name)) {
        return(object@enrichment)
    }
    object@enrichment[[name]]
})


#' @rdname qseaSlots
#' @export
setGeneric("setEnrichment", function(object, ...) {
    standardGeneric("setEnrichment")
})

#' @rdname qseaSlots
#' @export
setMethod("setEnrichment", "qseaSet", function(object, value, name = NULL) {
    if (is.null(name)) {
        object@enrichment <- value
    } else {
        object@enrichment[[name]] <- value
    }
    object
})


#' @rdname qseaSlots
#' @export
setGeneric("setCounts", function(object, ...) standardGeneric("setCounts"))

#' @rdname qseaSlots
#' @export
setMethod("setCounts", "qseaSet", function(object, value) {
    object@count_matrix <- as.matrix(value)
    object
})


#' @rdname qseaSlots
#' @export
setGeneric("setRegions", function(object, ...) standardGeneric("setRegions"))

#' @rdname qseaSlots
#' @export
setMethod("setRegions", "qseaSet", function(object, value) {
    object@regions <- value
    object
})


#' @rdname qseaSlots
#' @export
setGeneric("setSampleTable", function(object, ...) {
    standardGeneric("setSampleTable")
})

#' @rdname qseaSlots
#' @export
setMethod("setSampleTable", "qseaSet", function(object, value) {
    object@sampleTable <- as.data.frame(value)
    object
})


#' @rdname qseaSlots
#' @export
setGeneric("setParameters", function(object, ...) {
    standardGeneric("setParameters")
})

#' @rdname qseaSlots
#' @export
setMethod("setParameters", "qseaSet", function(object, values) {
    if (!is.list(values) || (length(values) > 0 && is.null(names(values)))) {
        stop("`values` must be a named list of parameters.")
    }
    object@parameters[names(values)] <- values
    object
})


#' @rdname qseaSlots
#' @export
setGeneric("setCNV", function(object, ...) standardGeneric("setCNV"))

#' @rdname qseaSlots
#' @export
setMethod("setCNV", "qseaSet", function(object, value) {
    object@cnv <- value
    object
})


# qsea::getCNV() subsets by sample name and returns NULL when the set has no
# CNV. Callers that need the slot exactly as stored, empty GRanges and all,
# use this instead.
.cnvSlot <- function(object) object@cnv


# Rewrite the sample names in every slot that carries them. The caller has
# already checked that `newNames` is valid, unique and in sample table order.
# The sample table goes last, because qsea::getCNV() looks the sample names up
# there and would not find the ones still stored in the other slots.
.setSampleNames <- function(object, newNames) {
    if (qsea::hasCNV(object)) {
        cnv <- qsea::getCNV(object)
        colnames(GenomicRanges::mcols(cnv)) <- newNames
        object <- setCNV(object, cnv)
    }

    zygosity <- qsea::getZygosity(object)
    rownames(zygosity) <- newNames
    object <- qsea::setZygosity(object, zygosity)

    counts <- qsea::getCounts(object)
    colnames(counts) <- newNames
    object <- setCounts(object, counts)

    for (type in c("file_name", "input_file")) {
        libraries <- getLibraries(object, type)
        if (!is.null(libraries)) {
            rownames(libraries) <- newNames
            object <- setLibraries(object, libraries, type)
        }
    }

    parameters <- getEnrichment(object, "parameters")
    if (!is.null(parameters)) {
        rownames(parameters) <- newNames
        object <- setEnrichment(object, parameters, "parameters")
    }

    factors <- getEnrichment(object, "factors")
    if (!is.null(factors)) {
        colnames(factors) <- newNames
        object <- setEnrichment(object, factors, "factors")
    }

    sampleTable <- qsea::getSampleTable(object) %>%
        tibble::remove_rownames() %>%
        dplyr::mutate(
            sample_name = newNames,
            rownameCol = newNames
        ) %>%
        tibble::column_to_rownames("rownameCol")

    setSampleTable(object, sampleTable)
}
