methods::setOldClass("prcomp")

# Sample (or sample group) IDs of one mesaPCA/mesaUMAP result: the row names
# of its coordinates.
.dimRedIds <- function(x) {
    if (methods::is(x, "mesaPCA")) {
        rownames(x@prcomp$x)
    } else {
        rownames(x@points)
    }
}

# Shared `windows` rules for mesaPCA and mesaUMAP; a message, or NULL if valid.
.checkWindows <- function(windows) {
    if (length(windows) == 0) return("`windows` must not be empty")
    if (anyNA(windows)) return("`windows` must not contain NA")
    if (anyDuplicated(windows)) return("`windows` must not contain duplicates")
    NULL
}

# ==============================
# mesaDimRed
# ==============================

#' Dimensionality reduction results container
#'
#' Aggregates one or more dimensionality reduction (DR) results (e.g., PCA/UMAP)
#' computed on a common set of samples.
#'
#' A valid object has `res` empty, or:
#' * `res` is a named list with unique names;
#' * its elements are all [mesaPCA-class] or all
#'   [mesaUMAP-class] objects, matching `params$method` if set;
#' * every element covers the same sample IDs, in the same order,
#'   with no duplicates;
#' * each sample ID is a value of `sampleTable$sample_name` (of
#'   `sampleTable$group` when `params$useGroupMeans` is `TRUE`);
#' * each sample ID is a column of `dataTable`, when it has columns.
#'
#' @section Accessors:
#' Read the slots with these functions rather than `@`:
#' * [getResults()]: the list of results (`res`);
#' * [getSampleTable()] and [getSampleNames()]: the sample annotations and the
#'   sample IDs the results cover;
#' * [getParameters()]: the parameters (`params`);
#' * [getDimRedData()]: the data used to compute the results (`dataTable`);
#' * [getCoordinates()] and [getWindowNames()]: per-result coordinates and
#'   window IDs, as named lists;
#' * [as.data.frame()][dimRedAccessors]: the coordinates of every result in
#'   one long table, with the sample annotations.
#'
#' @slot res `list`
#'   Individual DR result objects (e.g., [mesaPCA-class], [mesaUMAP-class]).
#'
#' The sample IDs of the results are the row names of each element of
#' `res` (`prcomp$x` for [mesaPCA-class],
#' `points` for [mesaUMAP-class]). With
#' `useGroupMeans = TRUE` they are sample group names.
#'
#' @slot sampleTable `data.frame`
#'   Sample annotations; the `sample_name` column holds the sample IDs
#'   (the row names mirror it, as in a `qseaSet`).
#'
#' @slot params `list`
#'   Parameters used to generate results in `res`.
#'
#' @slot dataTable `data.frame`
#'   Optional matrix/data used to compute DR (for reproducibility).
#'
#' @name mesaDimRed-class
#' @aliases mesaDimRed-class
#' @rdname mesaDimRed-class
#' @seealso [mesaPCA-class], [mesaUMAP-class], [dimRedAccessors],
#'   [getPCA()], [getUMAP()]
#' @exportClass mesaDimRed
setClass("mesaDimRed",
    slots = c(
        res = "list",
        sampleTable = "data.frame",
        params = "list",
        dataTable = "data.frame"
    )
)


#' Construct a `mesaDimRed` object
#'
#' @param res `list`
#'   DR result objects. **Default:** none. `res` must be supplied.
#'
#' @param sampleTable `data.frame`
#'   Sample annotations, with the sample IDs in a `sample_name` column.
#'   **Default:** none.
#'
#' @param params `list`
#'   Parameters used to compute `res`. **Default:** none.
#'
#' @param dataTable `data.frame`
#'   Optional numeric matrix/data used for DR.
#'   **Default:** `data.frame()`.
#'
#' @return A [mesaDimRed-class] object:
#' * stores DR results in `res`,
#' * carries sample metadata in `sampleTable`,
#' * and persists parameters/data in `params` / `dataTable`.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(20), nrow = 5, ncol = 4,
#'     dimnames = list(paste0("S", 1:5), paste0("W", 1:4)))
#' st <- data.frame(sample_name = rownames(x),
#'     group = rep(c("A", "B"), c(3, 2)), row.names = rownames(x))
#' mp <- mesaPCA(prcomp = stats::prcomp(x), windows = colnames(x))
#' md <- mesaDimRed(res = list(pca1 = mp), sampleTable = st,
#'     params = list(method = "PCA"))
#' md
#' getSampleNames(md)
#'
#' @rdname mesaDimRed-class
#' @export
mesaDimRed <- function(res, sampleTable, params, dataTable = data.frame()) {
    methods::new(
        "mesaDimRed",
        res = res, sampleTable = sampleTable,
        params = params, dataTable = dataTable
    )
}


#' @rdname mesaDimRed-class
#' @param object `mesaDimRed`
setMethod("show", "mesaDimRed", function(object) {
    cat("Object containing ", length(object@res),
        " dimensionality reduction objects for ",
        length(getSampleNames(object)), " samples", sep = "")
    cat("\n")
})

# Checks the sample IDs shared by every element of `res` against the other
# slots; a message, or NULL if valid.
.checkDimRedIds <- function(object, ids) {
    if (anyDuplicated(ids)) return("sample IDs in `res` must be unique")

    if (isTRUE(object@params$useGroupMeans)) {
        known <- object@sampleTable$group
    } else if ("sample_name" %in% colnames(object@sampleTable)) {
        known <- object@sampleTable$sample_name
    } else {
        return("`sampleTable` must have a `sample_name` column")
    }
    missing <- setdiff(ids, known)
    if (length(missing) > 0) {
        return(paste0("sample IDs in `res` are missing from `sampleTable`: ",
            paste(missing, collapse = ", ")))
    }

    if (ncol(object@dataTable) > 0) {
        missing <- setdiff(ids, colnames(object@dataTable))
        if (length(missing) > 0) {
            return(paste0("`dataTable` has no column for samples: ",
                paste(missing, collapse = ", ")))
        }
    }
    NULL
}

# Validity check
setValidity("mesaDimRed", function(object) {
    res <- object@res
    if (length(res) == 0) return(TRUE)

    isPCA <- vapply(res, methods::is, logical(1), "mesaPCA")
    isUMAP <- vapply(res, methods::is, logical(1), "mesaUMAP")
    if (!all(isPCA | isUMAP)) {
        return("every element of `res` must be a mesaPCA or mesaUMAP object")
    }
    if (!all(isPCA) && !all(isUMAP)) {
        return("`res` must not mix mesaPCA and mesaUMAP objects")
    }

    resNames <- names(res)
    if (is.null(resNames) || any(resNames == "") || anyDuplicated(resNames)) {
        return("`res` must be a named list with unique names")
    }

    method <- if (all(isPCA)) "PCA" else "UMAP"
    if (!is.null(object@params$method) &&
        !identical(object@params$method, method)) {
        return(paste0("`params$method` is \"", object@params$method,
            "\" but `res` holds ", method, " results"))
    }

    ids <- .dimRedIds(res[[1]])
    sameIds <- vapply(res, function(x) identical(.dimRedIds(x), ids),
        logical(1))
    if (!all(sameIds)) {
        return(paste("every element of `res` must cover the same samples,",
            "in the same order"))
    }
    msg <- .checkDimRedIds(object, ids)
    if (!is.null(msg)) return(msg)
    TRUE
})

setMethod("plotPCA", "mesaDimRed", plotPCA.mesaDimRed)

#' @rdname dimRedAccessors
#' @export
setMethod("getSampleTable", "mesaDimRed", function(object) {
    object@sampleTable
})

#' @rdname dimRedAccessors
setMethod("getSampleNames", "mesaDimRed", function(object) {
    if (length(object@res) == 0) return(character())
    .dimRedIds(object@res[[1]])
})

# ==============================
# mesaPCA
# ==============================

#' PCA results container
#'
#' Stores a PCA fit (from [stats::prcomp()]) computed over methylation windows.
#'
#' A valid object has non-empty `windows` with no `NA` or duplicates,
#' one per row of `prcomp$rotation`, and a `prcomp$x` matrix whose row
#' names are the sample IDs.
#'
#' @section Accessors:
#' * [getPrcomp()]: the `prcomp` fit, e.g. for `sdev` or `rotation`;
#' * [getCoordinates()]: the sample coordinates (`prcomp$x`), as a
#'   `data.frame`;
#' * [getWindowNames()]: the window IDs (`windows`).
#'
#' @slot prcomp `prcomp`
#'   A PCA fit returned by [stats::prcomp()].
#'
#' @slot windows `character()`
#'   Window IDs used in the PCA.
#'
#' @name mesaPCA-class
#' @aliases mesaPCA-class
#' @rdname mesaPCA-class
#' @seealso [mesaDimRed-class], [dimRedAccessors]
#' @exportClass mesaPCA
setClass("mesaPCA",
    slots = c(
        prcomp = "prcomp",
        windows = "character"
    )
)


#' Construct a `mesaPCA` object
#'
#' @param prcomp `prcomp`
#'   PCA fit from [stats::prcomp()]. **Default:** none.
#'
#' @param windows `character()`
#'   Window IDs used. **Default:** none.
#'
#' @return A [mesaPCA-class] object containing the PCA fit and window IDs.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(20), nrow = 5, ncol = 4,
#'     dimnames = list(paste0("S", 1:5), paste0("W", 1:4)))
#' pc <- stats::prcomp(x, center = TRUE, scale. = FALSE)
#' mp <- mesaPCA(prcomp = pc, windows = colnames(x))
#' mp
#'
#' @rdname mesaPCA-class
#' @export
mesaPCA <- function(prcomp, windows) {
    methods::new("mesaPCA", prcomp = prcomp, windows = windows)
}


#' @rdname mesaPCA-class
#' @param object `mesaPCA`
setMethod("show", "mesaPCA", function(object) {
    n <- tryCatch(nrow(object@prcomp$x), error = function(e) NA_integer_)
    cat("PCA result for ", n,
        " samples calculated over ", length(object@windows),
        " windows", sep = "")
    cat("\n")
})

# Validity check
setValidity("mesaPCA", function(object) {
    msg <- .checkWindows(object@windows)
    if (!is.null(msg)) return(msg)

    x <- object@prcomp$x
    if (!is.matrix(x) || is.null(rownames(x))) {
        return("`prcomp$x` must be a matrix with sample IDs as row names")
    }
    if (length(object@windows) != NROW(object@prcomp$rotation)) {
        return("`windows` must have one entry per row of `prcomp$rotation`")
    }
    TRUE
})



# ==============================
# mesaUMAP
# ==============================

#' UMAP results container
#'
#' Stores per-sample coordinates from a UMAP embedding computed over methylation
#' windows.
#'
#' A valid object has non-empty `windows` with no `NA` or duplicates,
#' and at least one row of numeric `points` with the sample IDs as row
#' names.
#'
#' @section Accessors:
#' * [getCoordinates()]: the sample coordinates (`points`);
#' * [getWindowNames()]: the window IDs (`windows`).
#'
#' @slot points `data.frame`
#'   One row per sample with UMAP coordinates (e.g., `UMAP1`, `UMAP2`).
#'   Row names are sample IDs.
#'
#' @slot windows `character()`
#'   Window IDs used to compute the embedding.
#'
#' @name mesaUMAP-class
#' @aliases mesaUMAP-class
#' @rdname mesaUMAP-class
#' @seealso [mesaDimRed-class], [dimRedAccessors]
#' @exportClass mesaUMAP

setClass("mesaUMAP",
    slots = c(
        points = "data.frame",
        windows = "character"
    )
)


#' Construct a `mesaUMAP` object
#'
#' @param points `data.frame`
#'   One row per sample with positions in UMAP space (row names = sample IDs).
#'   **Default:** none.
#'
#' @param windows `character()`
#'   Window IDs used. **Default:** none.
#'
#' @return A [mesaUMAP-class] object containing UMAP coordinates and window IDs.
#'
#' @examples
#' pts <- data.frame(UMAP1 = c(0.1, -0.2, 0.0),
#'     UMAP2 = c(0.3, 0.1, -0.1))
#' rownames(pts) <- paste0("S", 1:3)
#' mu <- mesaUMAP(points = pts, windows = c("w1", "w2", "w3"))
#' mu
#'
#' @rdname mesaUMAP-class
#' @export
mesaUMAP <- function(points, windows) {
    methods::new("mesaUMAP", points = points, windows = windows)
}


#' @rdname mesaUMAP-class
#' @param object `mesaUMAP`
setMethod("show", "mesaUMAP", function(object) {
    cat("UMAP result for ", nrow(object@points),
        " samples calculated over ", length(object@windows),
        " windows", sep = "")
    cat("\n")
})

# Validity check
setValidity("mesaUMAP", function(object) {
    msg <- .checkWindows(object@windows)
    if (!is.null(msg)) return(msg)

    if (nrow(object@points) == 0) return("`points` must have at least one row")
    if (!tibble::has_rownames(object@points)) {
        return("`points` must have sample IDs as row names")
    }
    if (!all(vapply(object@points, is.numeric, logical(1)))) {
        return("`points` columns must all be numeric")
    }
    TRUE
})


# ==============================
# Accessors
# ==============================

#' Accessors for dimensionality reduction results
#'
#' Read the contents of [mesaDimRed-class], [mesaPCA-class] and
#' [mesaUMAP-class] objects, as returned by [getPCA()] and [getUMAP()].
#'
#' @param object A [mesaDimRed-class] object, or a [mesaPCA-class] or
#'   [mesaUMAP-class] object where stated.
#'
#' @param x A [mesaDimRed-class] object.
#'
#' @param row.names,optional Not used.
#'
#' @param ... Not used.
#'
#' @return
#' * `getResults()`: a named list of [mesaPCA-class] or [mesaUMAP-class]
#'   objects, one per result (e.g. per `topVarNum` value); empty when there
#'   are no results.
#' * `getSampleTable()`: the sample annotations, a `data.frame` with the
#'   sample IDs in its `sample_name` column.
#' * `getSampleNames()`: the sample IDs the results cover (the row names of
#'   each result), or `character()` when there are no results. With
#'   `useGroupMeans = TRUE` they are sample group names.
#' * `getParameters()`: the list of parameters used to compute the results,
#'   including `method` (`"PCA"` or `"UMAP"`).
#' * `getDimRedData()`: the data used to compute the results (a
#'   `data.frame` of windows by samples), or an empty `data.frame` unless
#'   the results were made with `returnDataTable = TRUE`.
#' * `getCoordinates()`: for a [mesaPCA-class], its principal component
#'   scores (`prcomp$x`); for a [mesaUMAP-class], its UMAP coordinates. Each
#'   is a `data.frame` with one row per sample and sample IDs as row names.
#'   For a [mesaDimRed-class], a named list of these, one per result.
#' * `getPrcomp()`: the [stats::prcomp()] fit of a [mesaPCA-class].
#' * `as.data.frame()`: one `data.frame` with a row per result and sample:
#'   `resName` (the result's name, as in `getResults()`), `sample_name`, the
#'   coordinate columns (`PC1`, `PC2`, ... or `UMAP1`, `UMAP2`), then the
#'   `sampleTable` columns. With `useGroupMeans = TRUE` the ID column is
#'   `group` and no sample annotations are added.
#'
#' @seealso [getWindowNames()] for the window IDs of each result,
#'   [mesaDimRed-class], [mesaPCA-class], [mesaUMAP-class], [getPCA()],
#'   [getUMAP()], [plotDimRed()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' pca <- getPCA(exampleTumourNormal, topVarNum = c(100, 500),
#'     verbose = FALSE)
#'
#' names(getResults(pca))
#' head(getSampleTable(pca))
#' getSampleNames(pca)
#' getParameters(pca)$method
#'
#' first <- getResults(pca)[[1]]
#' head(getCoordinates(first))
#' summary(getPrcomp(first))$importance[, 1:3]
#'
#' # Coordinates of every result at once
#' lapply(getCoordinates(pca), dim)
#'
#' # Everything in one table, e.g. for a custom plot
#' head(as.data.frame(pca)[, 1:6])
#'
#' @name dimRedAccessors
NULL

#' @rdname dimRedAccessors
#' @export
setGeneric("getResults", function(object, ...) {
    standardGeneric("getResults")
})

#' @rdname dimRedAccessors
#' @export
setMethod("getResults", "mesaDimRed", function(object) object@res)

#' @rdname dimRedAccessors
#' @export
setMethod("getParameters", "mesaDimRed", function(object) object@params)

#' @rdname dimRedAccessors
#' @export
setGeneric("getDimRedData", function(object, ...) {
    standardGeneric("getDimRedData")
})

#' @rdname dimRedAccessors
#' @export
setMethod("getDimRedData", "mesaDimRed", function(object) object@dataTable)

#' @rdname dimRedAccessors
#' @export
setGeneric("getCoordinates", function(object, ...) {
    standardGeneric("getCoordinates")
})

#' @rdname dimRedAccessors
#' @export
setMethod("getCoordinates", "mesaPCA", function(object) {
    as.data.frame(object@prcomp$x)
})

#' @rdname dimRedAccessors
#' @export
setMethod("getCoordinates", "mesaUMAP", function(object) object@points)

#' @rdname dimRedAccessors
#' @export
setMethod("getCoordinates", "mesaDimRed", function(object) {
    lapply(object@res, getCoordinates)
})

#' @rdname dimRedAccessors
#' @export
setGeneric("getPrcomp", function(object, ...) {
    standardGeneric("getPrcomp")
})

#' @rdname dimRedAccessors
#' @export
setMethod("getPrcomp", "mesaPCA", function(object) object@prcomp)

#' @rdname dimRedAccessors
#' @importFrom BiocGenerics as.data.frame
#' @export
setMethod("as.data.frame", "mesaDimRed",
    function(x, row.names = NULL, optional = FALSE, ...) {
        byGroup <- isTRUE(getParameters(x)$useGroupMeans)
        idCol <- if (byGroup) "group" else "sample_name"
        if (length(getResults(x)) == 0) {
            out <- data.frame(resName = character(), id = character())
            return(rlang::set_names(out, c("resName", idCol)))
        }

        out <- getCoordinates(x) %>%
            lapply(tibble::rownames_to_column, var = idCol) %>%
            dplyr::bind_rows(.id = "resName")
        if (!byGroup) {
            # Coordinate columns keep their names if sampleTable reuses one.
            out <- dplyr::left_join(out, getSampleTable(x), by = idCol,
                suffix = c("", ".sampleTable"))
        }
        as.data.frame(out)
    }
)



# ==============================
# S3 methods for dplyr verbs
# ==============================

#' Mutate the sample table of a mesaDimRed
#'
#' Extend [dplyr::mutate()] to operate on the `sampleTable` slot of a
#' [mesaDimRed-class] object.
#'
#' @param .data `mesaDimRed`
#'   Object to modify.
#'
#' @param ...
#'   Arguments passed to [dplyr::mutate()].
#'
#' @return A `mesaDimRed` object:
#' * **sampleTable** updated with the mutated columns.
#' * Sample identity is preserved; attempts to alter `sample_name` are rejected.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' md <- exampleTumourNormal %>%
#'     getPCA() %>%
#'     mutate(group2 = paste0(group, "_2"))
#'
#' stopifnot("group2" %in% colnames(getSampleTable(md)))
#'
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom glue glue
#' @method mutate mesaDimRed
#' @export
mutate.mesaDimRed <- function(.data, ...) {
    newTable <- .data@sampleTable %>%
        tibble::rownames_to_column(".rownameCol") %>%
        dplyr::mutate(...) %>%
        tibble::column_to_rownames(".rownameCol")

    if (!identical(.data@sampleTable$sample_name, newTable$sample_name)) {
        stop(glue::glue(
            "sample_names cannot be changed with dplyr::mutate()."
        ))
    }

    .data@sampleTable <- newTable

    return(.data)}


#' Left-join onto the sample table of a mesaDimRed
#'
#' Extend [dplyr::left_join()] to operate on the `sampleTable` slot of a
#' [mesaDimRed-class] object.
#'
#' @param x `mesaDimRed`
#'   Object whose `sampleTable` will be joined.
#'
#' @param y `data.frame`
#'   Table to join.
#'
#' @param by `character()` or `NULL`
#'   Join columns; see [dplyr::left_join()]. **Default:** `NULL`.
#'
#' @param copy `logical(1)`
#'   See [dplyr::left_join()]. **Default:** `FALSE`.
#'
#' @param suffix `character(2)`
#'   Suffixes appended to overlapping non-join column names.
#'   **Default:** `c(".x", ".y")`.
#'
#' @param keep `logical(1)` or `NULL`
#'   Retain join keys from both tables. **Default:** `NULL`.
#'
#' @param ...
#'   Additional arguments to [dplyr::left_join()].
#'
#' @return A `mesaDimRed` object:
#' * **sampleTable** updated with `y` via a left join.
#'
#' @examples
#' new <- data.frame(sample_name = c("Colon1_N", "Colon2_N"), foo = c(1, 2))
#'
#' exampleTumourNormal %>%
#'     getPCA() %>%
#'     left_join(new, by = join_by(sample_name)) %>%
#'     getSampleTable()
#'
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @method left_join mesaDimRed
#' @export
left_join.mesaDimRed <- function(x, y, by = NULL, copy = FALSE,
                                    suffix = c(".x", ".y"), keep = NULL, ...) {
    x@sampleTable <- x@sampleTable %>%
        tibble::rownames_to_column("rownameCol") %>%
        dplyr::left_join(
            y, by = by, copy = copy, suffix = suffix, keep = keep, ...
        ) %>%
        tibble::column_to_rownames("rownameCol")

    return(x)}
