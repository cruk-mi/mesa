methods::setOldClass("prcomp")

#' @title The MesaDimRed class
#' @description
#' This is a class for containing the results of a dimensionality reduction on a qseaSet
#'
#' @slot res The results, a list containing individual PCA/UMAP objects.
#' @slot sampleTable A sampleTable from the qseaSet.
#' @slot samples A character vector of the samples used.
#' @slot params The parameters used during construction.
#' @slot dataTable The data used to construct the PCA/UMAP object (optional).
#' @export
setClass("mesaDimRed",
         slots = c(
           res = "list",
           sampleTable = "data.frame",
           samples = "character",
           params = "list",
           dataTable = "data.frame"
         )
)
#' @title Constructor function of a mesaDimRed object.
#' @description Constructor function of a mesaDimRed object.
#' @param res The results, a list containing individual PCA/UMAP objects.
#' @param sampleTable A sampleTable from the qseaSet.
#' @param samples A character vector of the samples used.
#' @param params The parameters used during construction.
#' @param dataTable The data used to construct the PCA/UMAP object (optional).
#' @return A `mesaDimRed` object.
#' @export
mesaDimRed <-  function(res, sampleTable, samples, params, dataTable = data.frame()) {
  methods::new("mesaDimRed", res = res, sampleTable = sampleTable, samples = samples, params = params, dataTable = dataTable)
}

#' @title The mesaPCA class
#' @description
#' This is a class for containing a PCA object calculated on a qseaSet
#'
#' @slot prcomp A principal component object output from [stats::prcomp()]
#' @slot windows A character vector of the windows used.
#' @export
setClass("mesaPCA",
         slots = c(
           prcomp = "prcomp",
           windows = "character"
         )
)
#' @title Constructor function of a mesaPCA object.
#' @description Constructor function of a mesaPCA object.
#' @param prcomp A principal component object output from [stats::prcomp()]
#' @param windows A character vector of the windows used.

#' @return A `mesaPCA` object.
#' @export
mesaPCA <-  function(prcomp, windows) {
  methods::new("mesaPCA", prcomp = prcomp, windows = windows)
}




#' @title The mesaUMAP class
#' @description
#' This is a class for containing UMAP objects calculated on a qseaSet
#'
#' @slot points Data frame containing one row for each sample, with the positions in UMAP space.
#' @slot windows A character vector of the windows used.
#' @export
setClass("mesaUMAP",
         slots = c(
           points = "data.frame",
           windows = "character"
         )
)
#' @title Constructor function of a mesaUMAP object.
#' @description Constructor function of a mesaUMAP object.
#' @param points Data frame containing one row for each sample, with the positions in UMAP space.
#' @param windows A character vector of the windows used.

#' @return A `mesaUMAP` object.
#' @export
mesaUMAP <-  function(points, windows) {
  methods::new("mesaUMAP", points = points, windows = windows)
}

setMethod("show", "mesaUMAP", function(object) {
  cat("UMAP result for ", nrow(object@points), "samples calculated over ", length(object@windows), " windows", sep = "")
})

setMethod("show", "mesaPCA", function(object) {
  cat("PCA result for ", nrow(object@prcomp$x), "samples calculated over ", length(object@windows), " windows", sep = "")
})

setMethod("show", "mesaDimRed", function(object) {
  cat("Object containing ", length(object@res)," dimensionality reduction objects for ", length(object@samples), " samples", sep = "")
})


setMethod("plotPCA", "mesaDimRed", plotPCA.mesaDimRed)

setMethod("getSampleTable", "mesaDimRed", function(object) {
  object@sampleTable
})

#' This function extends the dplyr function mutate to act on a mesaDimRed sampleTable.
#' @method mutate mesaDimRed
#' @importFrom dplyr mutate
#' @param .data A mesaDimRed to mutate
#' @param ... Other arguments to pass to dplyr::mutate
#' @return A mesaDimRed object with the sampleTable changed by a call to dplyr::mutate
#' @export
mutate.mesaDimRed <- function(.data, ...){

  newTable <- .data@sampleTable  %>%
    tibble::rownames_to_column(".rownameCol") %>%
    dplyr::mutate(...) %>%
    tibble::column_to_rownames(".rownameCol")

  if (!(identical(.data@sampleTable$sample_name, newTable$sample_name))) {
    stop(glue::glue("Error: sample_names cannot be changed with dplyr::mutate()."))
  }

  .data@sampleTable <- newTable

  return(.data)}

#' This function extends the dplyr function left_join to act on mesaDimRed sampleTable.
#' @method left_join mesaDimRed
#' @importFrom dplyr left_join
#' @param x A mesaDimRed to join data onto the sampleTable of
#' @param y A data frame to join with the sampleTable
#' @param by A character vector of variables to join by. If NULL, will perform a join on all common variables.
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x.
#' This allows you to join tables across srcs, but it is a potentially expensive operation so you must opt into it.
#' @param suffix 	If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2.
#' @param ... Other arguments to pass to dplyr::left_join
#' @param keep Should the join keys from both x and y be preserved in the output?
#' @return A mesaDimRed object with the sampleTable changed by a call to dplyr::left_join
#' @export
left_join.mesaDimRed <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x",".y"), keep = NULL, ...){

  x@sampleTable <- x@sampleTable %>%
    tibble::rownames_to_column("rownameCol") %>%
    dplyr::left_join(y, by = by, copy = copy, suffix = suffix, keep = keep,...) %>%
    tibble::column_to_rownames("rownameCol")

  return(x)}
