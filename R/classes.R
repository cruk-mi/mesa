methods::setOldClass("prcomp")

# ==============================
# mesaDimRed
# ==============================

#' Dimensionality reduction results container
#'
#' Aggregates one or more dimensionality reduction (DR) results (e.g., PCA/UMAP)
#' computed on a common set of samples.
#'
#' @slot res `list`  
#'   Individual DR result objects (e.g., [mesaPCA-class], [mesaUMAP-class]).
#'
#' @slot sampleTable `data.frame`  
#'   Sample annotations; row names are sample IDs.
#'
#' @slot samples `character()`  
#'   Vector of sample IDs used.
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
#' @seealso [mesaPCA-class], [mesaUMAP-class]
#' @exportClass mesaDimRed
setClass("mesaDimRed",
         slots = c(
           res = "list",
           sampleTable = "data.frame",
           samples = "character",
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
#'   Sample annotations; row names must be sample IDs. **Default:** none.
#'
#' @param samples `character()`  
#'   Sample IDs included. **Default:** none.
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
#' * records the samples in `samples`,
#' * and persists parameters/data in `params` / `dataTable`.
#'
#' @examples
#' st <- data.frame(sample_name = c("S1","S2"),
#'                  group = c("A","B"),
#'                  row.names = "sample_name")
#' md <- mesaDimRed(res = list(), sampleTable = st,
#'                  samples = rownames(st), params = list(),
#'                  dataTable = data.frame())
#' md
#'
#' @rdname mesaDimRed-class
#' @export
mesaDimRed <-  function(res, sampleTable, samples, params, dataTable = data.frame()) {
  methods::new("mesaDimRed", res = res, sampleTable = sampleTable, samples = samples, params = params, dataTable = dataTable)
}


#' @rdname mesaDimRed-class
#' @param object `mesaDimRed`  
setMethod("show", "mesaDimRed", function(object) {
  cat("Object containing ", length(object@res),
      " dimensionality reduction objects for ",
      length(object@samples), " samples", sep = "")
  cat("\n")
})

# Validity check
setValidity("mesaDimRed", function(object) {
  if (!is.data.frame(object@sampleTable)) return("`sampleTable` must be a data.frame")
  if (!is.character(object@samples)) return("`samples` must be character")
  TRUE
})

setMethod("plotPCA", "mesaDimRed", plotPCA.mesaDimRed)

setMethod("getSampleTable", "mesaDimRed", function(object) {
  object@sampleTable
})

# ==============================
# mesaPCA
# ==============================

#' PCA results container
#'
#' Stores a PCA fit (from [stats::prcomp()]) computed over methylation windows.
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
#' @seealso [mesaDimRed-class]
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
#'             dimnames = list(paste0("S",1:5), paste0("W",1:4)))
#' pc <- stats::prcomp(x, center = TRUE, scale. = FALSE)
#' mp <- mesaPCA(prcomp = pc, windows = colnames(x))
#' mp
#'
#' @rdname mesaPCA-class
#' @export
mesaPCA <-  function(prcomp, windows) {
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
  if (!inherits(object@prcomp, "prcomp")) return("`prcomp` must be a stats::prcomp object")
  if (!is.character(object@windows)) return("`windows` must be character")
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
#' @seealso [mesaDimRed-class]
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
#'                   UMAP2 = c(0.3,  0.1, -0.1))
#' rownames(pts) <- paste0("S", 1:3)
#' mu <- mesaUMAP(points = pts, windows = c("w1","w2","w3"))
#' mu
#'
#' @rdname mesaUMAP-class
#' @export
mesaUMAP <-  function(points, windows) {
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
  if (!is.data.frame(object@points)) return("`points` must be a data.frame")
  if (!is.character(object@windows)) return("`windows` must be character")
  TRUE
})



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
#' st <- data.frame(sample_name = c("A","B"),
#'                  group = c("X","Y"),
#'                  row.names = "sample_name")
#' md <- mesaDimRed(res = list(), sampleTable = st,
#'                  samples = rownames(st), params = list())
#' md2 <- mutate(md, group2 = paste0(group, "_2"))
#' stopifnot("group2" %in% colnames(
#'   methods::slot(md2, "sampleTable")
#' ))
#'
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom glue glue
#' @method mutate mesaDimRed
#' @export
mutate.mesaDimRed <- function(.data, ...) {
  newTable <- .data@sampleTable  %>%
    tibble::rownames_to_column(".rownameCol") %>%
    dplyr::mutate(...) %>%
    tibble::column_to_rownames(".rownameCol")
  
  if (!identical(.data@sampleTable$sample_name, newTable$sample_name)) {
    stop(glue::glue("Error: sample_names cannot be changed with dplyr::mutate()."))
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
#' st <- data.frame(sample_name = c("A","B"),
#'                  group = c("X","Y"),
#'                  row.names = "sample_name")
#' md <- mesaDimRed(res = list(), sampleTable = st,
#'                  samples = rownames(st), params = list())
#' ann <- data.frame(rownameCol = c("A","B"), batch = c(1,2))
#' md3 <- left_join(md, ann, by = "rownameCol")
#'
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @method left_join mesaDimRed
#' @export
left_join.mesaDimRed <- function(x, y, by = NULL, copy = FALSE,
                                 suffix = c(".x",".y"), keep = NULL, ...) {
  x@sampleTable <- x@sampleTable %>%
    tibble::rownames_to_column("rownameCol") %>%
    dplyr::left_join(y, by = by, copy = copy, suffix = suffix, keep = keep, ...) %>%
    tibble::column_to_rownames("rownameCol")
  
  return(x)}
