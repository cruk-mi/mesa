#' This function extends the dplyr function filter to act on qseaSet sampleTable.
#' @method filter qseaSet
#' @importFrom dplyr filter
#' @param qseaSet A qseaSet to filter on
#' @param ... Other arguments to pass to dplyr::filter
#' @param .preserve Not implemented as irrelevant for qseaSet
#' @return A qseaSet object with the sampleTable enhanced with the information on number of reads etc
#' @export filter.qseaSet
#' @export
filter.qseaSet <- function(.data, ..., .preserve = FALSE){filterQset(.data, ...)}

#' This function takes a qseaSet and filters the samples in it based on a call to dplyr::filter on the sampleTable.
#' @param qseaSet The qseaSet object.
#' @param ... Other arguments to pass to dplyr::filter
#' @return A qseaSet object with the sampleTable enhanced with the information on number of reads etc
#' @export
filterQset <- function(qseaSet, ...){
  namesToKeep <- qseaSet %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(...) %>%
    dplyr::pull(sample_name)

  if ("PooledControl" %in% qsea::getSampleNames(qseaSet)) {
    namesToKeep <- union(namesToKeep, "PooledControl")
  }

  subsetQset(qseaSet, samplesToKeep = namesToKeep) %>%
    return()
}

#' This function extends the dplyr function mutate to act on qseaSet sampleTable.
#' @method mutate qseaSet
#' @importFrom dplyr mutate
#' @param .data A qseaSet to mutate
#' @param ... Other arguments to pass to dplyr::mutate
#' @return A qseaSet object with the sampleTable changed by a call to dplyr::mutate
#' @export mutate.qseaSet
#' @export
mutate.qseaSet <- function(.data, ...){mutateQset(.data, ...)}


#' This function takes a qseaSet and mutates its sampleTable based on a call to dplyr::mutate.
#' @param qseaSet The qseaSet object.
#' @param ... Other arguments to pass to dplyr::mutate
#' @return A qseaSet object with the sampleTable mutated as specified.
#' @export
mutateQset <- function(qseaSet, ...){
  qseaSet@sampleTable <- qseaSet %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column("rownameCol") %>%
    dplyr::mutate(...) %>%
    tibble::column_to_rownames("rownameCol")

  return(qseaSet)
}

#' This function extends the dplyr function left_join to act on qseaSet sampleTable.
#' @method left_join qseaSet
#' @importFrom dplyr left_join
#' @param qseaSet A qseaSet to mutate
#' @param ... Other arguments to pass to dplyr::mutate
#' @return A qseaSet object with the sampleTable changed by a call to dplyr::left_join
#' @export
left_join.qseaSet <- function(qseaSet, ...){leftJoinQset(qseaSet, ...)}

#' This function takes a qseaSet and appends data to its sampleTable based on a call to dplyr::left_join
#' @param qseaSet The qseaSet object.
#' @param newData The new data frame to join the sampleTable on.
#' @param ... Other arguments to pass to dplyr::left_join. Can include "by" to specify which columns to join on
#' @return A qseaSet object with the sampleTable appended to as specified.
#' @export
leftJoinQset <- function(qseaSet, newData, ...){
  qseaSet@sampleTable <- qseaSet %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column("rownameCol") %>%
    dplyr::left_join(newData,...) %>%
    tibble::column_to_rownames("rownameCol")

  return(qseaSet)
}

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
  qseaSet@sampleTable <- qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::select(...)

  return(qseaSet)
}

#' This function extends the dplyr function pull to act on qseaSet sampleTable.
#' @method pull qseaSet
#' @importFrom dplyr pull
#' @param qseaSet A qseaSet to pull a column from
#' @param var A variable specified as a column name or an integer (negative counting from right)
#' @param name An optional parameter that specifies the column to be used as names for a named vector. Specified in a similar manner as var.
#' @param ... Other arguments to pass to dplyr::pull (column name to extract)
#' @return A vector with the contents of the column of the sample table that was extracted.
#' @export pull.qseaSet
#' @export
pull.qseaSet <- function(.data, ...){pullQset(.data, ...)}


#' This function takes a qseaSet and pulls a column from its sampleTable based on a call to dplyr::pull
#' @param qseaSet The qseaSet object.
#' @param ... Other arguments to pass to dplyr::pull
#' @return A column of the sampleTable as a vector
#' @export
pullQset <- function(qseaSet, ...){
  qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::pull(...) %>%
    return()
}

#' This function extends the function sort to act on qseaSet, to reorder the names. It uses gtools::mixedsort to sort numbers correctly (i.e. 10 does not come before 2).
#' @method sort qseaSet
#' @param x A qseaSet to reorder the samples of.
#' @param decreasing Whether to order the qseaSet in reverse alphabetical order
#' @param ... Other arguments to pass to gtools::mixedsort
#' @return A qseaSet with the sample names reordered alphabetically
#' @export sort.qseaSet
#' @export
sort.qseaSet <- function(x, ...){sortQset(x, decreasing = FALSE, ...)}


#' This function takes a qseaSet and sorts the samples in it, so they are in alphabetical order (with numbers increasing via mixedsort)
#' @param qseaSet The qseaSet object.
#' @param decreasing Whether to order the qseaSet in reverse alphabetical order
#' @return A qseaSet object with the objects sorted
#' @export
sortQset <- function(qseaSet, decreasing = FALSE){
  qseaSet %>%
    subsetQset(samplesToKeep = (qseaSet %>% qsea::getSampleNames() %>% gtools::mixedsort(decreasing = decreasing))) %>%
    return()

}

#' This function takes a qseaSet object and filters the regions inside it by a call to dplyr::filter.
#' @param qseaSet The qseaSet object.
#' @param ... Additional arguments to be used to filter the regions ONLY, as if they were a data frame.
#' @return A qseaSet object, with the regions filtered appropriately.
#' @export
filterRegions <- function(qseaSet, ...){

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
#' @param qseaSet A qseaSet to reorder the samples in.
#' @param ... Other arguments to pass to dplyr::arrange
#' @param .by_group Not implemented as qsea requires a data frame not a tibble for the sampleTable.
#' @return A qseaSet with the sample names reordered according to a column of the sampleTable.
#' @export arrange.qseaSet
#' @export
arrange.qseaSet <- function(.data, ..., .by_group = FALSE){arrangeQset(.data, ..., .by_group)}

#' This function takes a qseaSet and sorts the samples inside it using dplyr::arrange, using a column(s) of the sample table. Can also use desc() to reverse sort, as in dplyr::arrange.
#' @param qseaSet The qseaSet object.
#' @param ... Other arguments to pass to dplyr::arrange
#' @param .by_group Not implemented as qsea requires a data frame not a tibble for the sampleTable.
#' @return A qseaSet object with the samples sorted according to the column(s) selected.
#' @export
arrangeQset <- function (qseaSet, ..., .by_group = FALSE) {

  sampleOrder <- qseaSet %>%
    qsea::getSampleTable() %>%
    dplyr::arrange(...) %>%
    dplyr::pull(sample_name)

  qseaSet %>%
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
