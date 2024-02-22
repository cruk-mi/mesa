#' This function extends the dplyr function filter to act on qseaSet, subseting samples by using the sampleTable.
#' @method filter qseaSet
#' @importFrom dplyr filter
#' @param .data A qseaSet to filter, based on the sampleTable
#' @param ... Other arguments to pass to dplyr::filter
#' @param .preserve Not implemented as grouping is not available  for qseaSet
#' @return A qseaSet object, with only samples that are selected by the filtering operation.
#' @export
filter.qseaSet <- function(.data, ..., .preserve = FALSE){

  namesToKeep <- .data %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(...) %>%
    dplyr::pull(sample_name)

  return(subsetQset(.data, samplesToKeep = namesToKeep))

  }

#' This function extends the dplyr function mutate to act on qseaSet sampleTable.
#' @method mutate qseaSet
#' @importFrom dplyr mutate
#' @param .data A qseaSet to mutate
#' @param ... Other arguments to pass to dplyr::mutate
#' @return A qseaSet object with the sampleTable changed by a call to dplyr::mutate
#' @export
mutate.qseaSet <- function(.data, ...){

   newTable <- .data %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column(".rownameCol") %>%
    dplyr::mutate(...) %>%
    tibble::column_to_rownames(".rownameCol")
  
  if (!(identical(.data@sampleTable$sample_name, newTable$sample_name))) {
    stop(glue::glue("Error: sample_name cannot be changed with dplyr::mutate(). Use mesa::renameQsetNames() or mesa::renameSamples() instead."))
  }

  .data@sampleTable <- newTable

  if (!(identical(.data@sampleTable$sample_name, newTable$sample_name))) {
    stop(glue::glue("Error: sample_name cannot be changed with dplyr::mutate(). Use mesa::renameQsetNames() or mesa::renameSamples() instead."))
  }

  .data@sampleTable <- newTable

  return(.data)}

#' This function extends the dplyr function left_join to act on qseaSet sampleTable.
#' @method left_join qseaSet
#' @importFrom dplyr left_join
#' @param x A qseaSet to join data onto the sampleTable of
#' @param y A data frame to join with the sampleTable
#' @param by A character vector of variables to join by. If NULL, will perform a join on all common variables.
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x.
#' This allows you to join tables across srcs, but it is a potentially expensive operation so you must opt into it.
#' @param suffix 	If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2.
#' @param ... Other arguments to pass to dplyr::left_join
#' @param keep Should the join keys from both x and y be preserved in the output?
#' @return A qseaSet object with the sampleTable changed by a call to dplyr::left_join
#' @export
left_join.qseaSet <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x",".y"), keep = NULL, ...){

  x@sampleTable <- x %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column("rownameCol") %>%
    dplyr::left_join(y, by = by, copy = copy, suffix = suffix, keep = keep,...) %>%
    tibble::column_to_rownames("rownameCol")

  return(x)}

#' This function extends the dplyr function select to act on qseaSet sampleTable. Can also be used to rename columns.
#' @method select qseaSet
#' @importFrom dplyr select
#' @param .data A qseaSet to select columns of the sampleTable from
#' @param ... Other arguments to pass to dplyr::select
#' @return A qseaSet object with the sampleTable changed by a call to dplyr::select
#' @export select.qseaSet
#' @export
select.qseaSet <- function(.data, ...){selectQset(.data, ...)}

#' This function takes a qseaSet and selects columns from its sampleTable based on a call to dplyr::select
#' @param qseaSet The qseaSet object.
#' @param ... Other arguments to pass to dplyr::select
#' @return A qseaSet object with the sampleTable selected as specified.
#' @export
selectQset <- function(qseaSet, ...){
  newSampleTable <- qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::select(...)

  #ensure that the sample_name and group columns are still present at the front of the output
  qseaSet@sampleTable <- qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::select(sample_name, group) %>%
    dplyr::bind_cols(newSampleTable %>% dplyr::select(-tidyselect::matches("^sample_name$|^group$")))

  return(qseaSet)
}

#' This function extends the dplyr function pull to act on qseaSet sampleTable. It returns a column of the sampleTable.
#' @method pull qseaSet
#' @importFrom dplyr pull
#' @param .data A qseaSet.
#' @param var A variable specified as a column name or an integer (negative counting from right)
#' @param name An optional parameter that specifies the column to be used as names for a named vector. Specified in a similar manner as var.
#' @param ... Other arguments to pass to dplyr::pull (column name to extract)
#' @return A vector the same size as the number of samples in the qseaSet.
#' @export pull.qseaSet
#' @export
pull.qseaSet <- function(.data, var = -1, name = NULL, ...){

  .data <- .data %>% qsea::getSampleTable()

  #code copied directly from dplyr::pull
  var <- tidyselect::vars_pull(names(.data), !!rlang::enquo(var))
  name <- rlang::enquo(name)
  if (rlang::quo_is_null(name)) {
    return(.data[[var]])
  }
  name <- tidyselect::vars_pull(names(.data), !!name)
  rlang::set_names(.data[[var]], nm = .data[[name]])

  }

#' This function extends the function sort to act on qseaSet, to reorder the names. It uses gtools::mixedsort to sort numbers correctly (i.e. 10 does not come before 2).
#' @method sort qseaSet
#' @param x A qseaSet to reorder the samples of.
#' @param decreasing Whether to order the qseaSet in reverse alphabetical order
#' @param ... Other arguments to pass to gtools::mixedsort
#' @return A qseaSet with the sample names reordered alphabetically
#' #' @examples
#' sort(exampleTumourNormal, decreasing = TRUE)
#' @export
sort.qseaSet <- function(x, decreasing = FALSE, ...){

  x %>%
    subsetQset(samplesToKeep = (x %>% qsea::getSampleNames() %>% gtools::mixedsort(decreasing = decreasing))) %>%
    return()

  }

#' This function takes a qseaSet object and filters the regions inside it by a call to dplyr::filter.
#' @param qseaSet The qseaSet object.
#' @param ... Additional arguments to be used to filter the regions ONLY, as if they were a data frame.
#' @return A qseaSet object, with the regions filtered appropriately.
#' @export
filterWindows <- function(qseaSet, ...){

  reducedRegions <- qseaSet %>%
    qsea::getRegions() %>%
    tibble::as_tibble() %>%
    dplyr::filter(...) %>%
    plyranges::as_granges()

  qseaSet %>%
    filterByOverlaps(reducedRegions) %>%
    return()

}

#' This function extends the function dplyr::arrange to act on qseaSet, to reorder the samples in the qseaSet.
#' @method arrange qseaSet
#' @importFrom dplyr arrange
#' @param .data A qseaSet to reorder the samples in.
#' @param ... Other arguments to pass to dplyr::arrange
#' @param .by_group Not implemented as qsea requires a data frame not a tibble for the sampleTable.
#' @return A qseaSet with the sample names reordered according to a column of the sampleTable.
#' @export arrange.qseaSet
#' @export
arrange.qseaSet <- function(.data, ..., .by_group = FALSE){

  sampleOrder <- .data %>%
    qsea::getSampleTable() %>%
    dplyr::arrange(...) %>%
    dplyr::pull(sample_name)

  .data %>%
    subsetQset(samplesToKeep = sampleOrder) %>%
    return()

  }

#' This function extends the base function length to return the number of samples in the qseaSet
#' @method colnames qseaSet
#' @param qseaSet A qseaSet to return the number of samples in
#' @return The number of samples in the qseaSet
#' @export colnames.qseaSet
#' @export
colnames.qseaSet <- function(qseaSet){colnames(qsea::getSampleTable(qseaSet))}
