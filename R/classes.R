methods::setOldClass("prcomp")

# ==============================
# mesaDimRed
# ==============================

#' Dimensionality reduction results container
#'
#' Aggregates one or more dimensionality reduction results (e.g., PCA/UMAP)
#' computed on a common set of samples.
#'
#' @slot res A list of individual DR results (e.g., \linkS4class{mesaPCA},
#'   \linkS4class{mesaUMAP}, or related objects).
#' @slot sampleTable A data.frame describing samples (row names = sample IDs).
#' @slot samples Character vector of sample IDs used.
#' @slot params List of parameters used to generate results in \code{res}.
#' @slot dataTable Optional data.frame containing the matrix used to compute DR.
#'
#' @name mesaDimRed-class
#' @aliases mesaDimRed-class
#' @rdname mesaDimRed-class
#' @seealso \linkS4class{mesaPCA}, \linkS4class{mesaUMAP}
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

#' Construct a \code{mesaDimRed} object
#'
#' @param res List of DR result objects.
#' @param sampleTable Data.frame with sample annotations (row names = sample IDs).
#' @param samples Character vector of sample IDs included.
#' @param params List of parameters used to compute \code{res}.
#' @param dataTable Optional data.frame of the numeric matrix used for DR.
#'
#' @return A \linkS4class{mesaDimRed} object.
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
#' @param object A \linkS4class{mesaDimRed} object.
#' @exportMethod show
setMethod("show", "mesaDimRed", function(object) {
  cat("Object containing ", length(object@res),
      " dimensionality reduction objects for ",
      length(object@samples), " samples", sep = "")
  cat("\n")
})

# (Optional) simple validity
setValidity("mesaDimRed", function(object) {
  if (!is.data.frame(object@sampleTable)) return("`sampleTable` must be a data.frame")
  if (!is.character(object@samples)) return("`samples` must be character")
  TRUE
})


# ==============================
# mesaPCA
# ==============================

#' PCA results container
#'
#' Stores a PCA fit (from \code{stats::prcomp()}) computed over methylation
#' windows.
#'
#' @slot prcomp A \code{stats::prcomp} object.
#' @slot windows Character vector of window IDs used in the PCA.
#'
#' @name mesaPCA-class
#' @aliases mesaPCA-class
#' @rdname mesaPCA-class
#' @seealso \linkS4class{mesaDimRed}
#' @exportClass mesaPCA
setClass("mesaPCA",
         slots = c(
           prcomp = "prcomp",
           windows = "character"
         )
)

#' Construct a \code{mesaPCA} object
#'
#' @param prcomp A PCA fit returned by \code{stats::prcomp()}.
#' @param windows Character vector of window IDs used.
#'
#' @return A \linkS4class{mesaPCA} object.
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
#' @param object A \linkS4class{mesaPCA} object.
#' @exportMethod show
setMethod("show", "mesaPCA", function(object) {
  n <- tryCatch(nrow(object@prcomp$x), error = function(e) NA_integer_)
  cat("PCA result for ", n,
      " samples calculated over ", length(object@windows),
      " windows", sep = "")
  cat("\n")
})

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
#' @slot points A data.frame with one row per sample and UMAP coordinates
#'   (e.g., columns `UMAP1`, `UMAP2`). Row names are sample IDs.
#' @slot windows Character vector of window IDs used to compute the embedding.
#'
#' @name mesaUMAP-class
#' @aliases mesaUMAP-class
#' @rdname mesaUMAP-class
#' @seealso \linkS4class{mesaDimRed}
#' @exportClass mesaUMAP
setClass("mesaUMAP",
         slots = c(
           points = "data.frame",
           windows = "character"
         )
)

#' Construct a \code{mesaUMAP} object
#'
#' @param points A data.frame with one row per sample with the positions in UMAP space. 
#'   Row names should be sample IDs.
#' @param windows Character vector of window IDs used.
#'
#' @return A \linkS4class{mesaUMAP} object.
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
#' @param object A \linkS4class{mesaUMAP} object.
#' @exportMethod show
setMethod("show", "mesaUMAP", function(object) {
  cat("UMAP result for ", nrow(object@points),
      " samples calculated over ", length(object@windows),
      " windows", sep = "")
  cat("\n")
})

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
#' Extends \code{dplyr::mutate()} to operate on the \code{sampleTable} slot.
#'
#' @param .data A \linkS4class{mesaDimRed} object.
#' @param ... Arguments passed to \code{dplyr::mutate()}.
#'
#' @return A \linkS4class{mesaDimRed} with an updated \code{sampleTable}.
#' @examples
#' st <- data.frame(sample_name = c("A","B"),
#'                  group = c("X","Y"),
#'                  row.names = "sample_name")
#' md <- mesaDimRed(res = list(), sampleTable = st,
#'                  samples = rownames(st), params = list())
#' md2 <- mutate(md, group2 = paste0(group, "_2"))
#' stopifnot("group2" %in% colnames(md2@sampleTable))
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
#' Extends \code{dplyr::left_join()} to operate on the \code{sampleTable} slot.
#'
#' @param x A \linkS4class{mesaDimRed}.
#' @param y A data.frame to join.
#' @param by,suffix,copy,keep,... Passed to \code{dplyr::left_join()}.
#'
#' @return A \linkS4class{mesaDimRed} with an updated \code{sampleTable}.
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
