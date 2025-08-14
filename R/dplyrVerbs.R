#' Filter samples in a qseaSet
#'
#' Extends [dplyr::filter()] to subset **samples** in a `qseaSet` using
#' predicates applied to its `sampleTable`.
#'
#' @details
#' The filtering is performed on `qsea::getSampleTable(.data)`, and the
#' selected `sample_name`s are used to subset the `qseaSet`. Grouped filtering
#' is not supported for `qseaSet`.
#'
#' @method filter qseaSet
#' @importFrom dplyr filter
#' @param .data A `qseaSet` to filter (by rows of its `sampleTable`).
#' @param ...  Predicate expressions passed to [dplyr::filter()].
#' @param .preserve Ignored; grouping is not implemented for `qseaSet`.
#'
#' @return The input qseaSet object, with only samples selected by the filtering operation.
#' @seealso [qsea::getSampleTable()], [subsetQset()], [selectQset()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   # keep only older patients
#'   exampleTumourNormal %>% filter(age >= 70)
#' }
#' }
#' @export
filter.qseaSet <- function(.data, ..., .preserve = FALSE){

  namesToKeep <- .data %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(...) %>%
    dplyr::pull(sample_name)

  return(subsetQset(.data, samplesToKeep = namesToKeep))

  }

#' Mutate columns in a qseaSet sample table
#'
#' Extends [dplyr::mutate()] to modify or add columns in a `qseaSet`'s
#' `sampleTable`.
#'
#' @details
#' The order of samples is preserved. The `sample_name` column **cannot be
#' changed** here; use `renameQsetNames()`/`renameSamples()` instead.
#'
#' @method mutate qseaSet
#' @importFrom dplyr mutate
#' @param .data A `qseaSet`.
#' @param ...  Mutations passed to [dplyr::mutate()].
#'
#' @return The input `qseaSet` object with an updated `sampleTable`.
#'
#' @seealso [qsea::getSampleTable()], [renameSamples()], [renameQsetNames()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   exampleTumourNormal %>% mutate(over70 = age > 70)
#' }
#' }
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

#' Left join data onto a qseaSet sample table
#'
#' Extends [dplyr::left_join()] to merge a data frame `y` into a `qseaSet`'s
#' `sampleTable`.
#'
#' @method left_join qseaSet
#' @importFrom dplyr left_join
#' @param x A `qseaSet` whose `sampleTable` will receive new columns.
#' @param y A data frame to join.
#' @param by Variables to join by (see [dplyr::left_join()]). If NULL, will perform a join on all common variables.
#' @param by Character or named character vector of join keys (see [dplyr::left_join()]).
#'   For `qseaSet`, this is typically `"sample_name"`. If omitted, dplyr uses the
#'   intersection of column names; specifying `by` explicitly is safer.
#' @param copy Logical; **remote backends only** (e.g., databases via dbplyr).
#'   When `TRUE`, `y` is copied into the same data source as `x` so the join can run.
#'   This can be expensive. Ignored for in-memory data frames/tibbles.
#' @param suffix Length-2 character vector with suffixes appended to non-joined
#'   duplicate column names from `x` and `y` (default `c(".x", ".y")`).
#' @param keep Logical; if `TRUE`, retain the join keys from both `x` and `y`
#'   in the output. If `FALSE`, keep only keys from `x` (dplyr default).
#' @param ... Additional arguments passed to [dplyr::left_join()], e.g.
#'   `relationship = "one-to-one"` to assert key uniqueness or
#'   `unmatched = "error"` to fail when keys don’t match (dplyr ≥ 1.1).
#'   
#' @return The input `qseaSet` with an updated `sampleTable`.
#'
#' @seealso [qsea::getSampleTable()], [dplyr::left_join()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   meta <- data.frame(sample_name = qsea::getSampleNames(exampleTumourNormal),
#'                      batch = rep(1:2, length.out = length(qsea::getSampleNames(exampleTumourNormal))))
#'   exampleTumourNormal %>% left_join(meta, by = "sample_name")
#' }
#' }
#' @export
left_join.qseaSet <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x",".y"), keep = NULL, ...){

  x@sampleTable <- x %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column("rownameCol") %>%
    dplyr::left_join(y, by = by, copy = copy, suffix = suffix, keep = keep,...) %>%
    tibble::column_to_rownames("rownameCol")

  return(x)}

#' Select/rename columns in a qseaSet sample table
#'
#' Extends [dplyr::select()] to keep or rename columns in a `qseaSet`'s
#' `sampleTable`.
#'
#' @details
#' This is a thin wrapper around [selectQset()], which preserves `sample_name`
#' and `group` at the front of the table if present.
#'
#' @method select qseaSet
#' @importFrom dplyr select
#' @param .data A `qseaSet`.
#' @param ...  Selection helpers passed to [dplyr::select()].
#'
#' @return The input `qseaSet` object with an updated `sampleTable`.
#'
#' @seealso [selectQset()], [qsea::getSampleTable()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   exampleTumourNormal %>% select(sample_name, group, age)
#' }
#' }
#' @export
select.qseaSet <- function(.data, ...){selectQset(.data, ...)}

#' Select columns of a qseaSet sample table
#'
#' Keeps/renames columns from the `sampleTable` of a `qseaSet`. If present,
#' `sample_name` and `group` are kept at the front of the table.
#'
#' @param qseaSet A `qseaSet`.
#' @param ...  Selection helpers passed to [dplyr::select()].
#'
#' @return The input `qseaSet` object with `sampleTable` updated.
#'
#' @seealso [select.qseaSet()], [qsea::getSampleTable()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   selectQset(exampleTumourNormal, sample_name, group)
#' }
#' }
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

#' Pull a column from a qseaSet sample table
#'
#' Extends [dplyr::pull()] to extract a column from the `sampleTable` of a `qseaSet`.
#'
#' @method pull qseaSet
#' @importFrom dplyr pull
#' @param .data A `qseaSet`.
#' @param var  Column to extract (name or position). Defaults to last column.
#' @param name Optional column to use for names in the returned vector.
#' @param ...  Passed to [dplyr::pull()].
#'
#' @return A vector with length equal to the number of samples in the qseaSet.
#'
#' @seealso [qsea::getSampleTable()], [dplyr::pull()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   exampleTumourNormal %>% pull(group)
#' }
#' }
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

#' Sort samples in a qseaSet
#'
#' Extends [base::sort()] to reorder the **samples** of a `qseaSet` using
#' [gtools::mixedsort()] so that numeric parts sort naturally (e.g., `"2"` before `"10"`).
#'
#' @method sort qseaSet
#' @param x A `qseaSet` to reorder.
#' @param decreasing Logical; reverse the order.
#' @param ... Passed to [gtools::mixedsort()].
#'
#' @return The input `qseaSet` object with samples reordered.
#'
#' @seealso [qsea::getSampleNames()], [subsetQset()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   sort(exampleTumourNormal, decreasing = TRUE)
#' }
#' }
#' @export
sort.qseaSet <- function(x, decreasing = FALSE, ...){

  x %>%
    subsetQset(samplesToKeep = (x %>% qsea::getSampleNames() %>% gtools::mixedsort(decreasing = decreasing))) %>%
    return()

  }

#' Filter regions (windows) inside a qseaSet
#'
#' Filters the **regions** of a `qseaSet` using [dplyr::filter()] predicates
#' applied to its regions as a data frame, then subsets the object accordingly.
#'
#' @param qseaSet A `qseaSet`.
#' @param ... Predicates passed to [dplyr::filter()] (applied to the regions).
#'
#' @return The input `qseaSet` object with regions filtered.
#'
#' @seealso [qsea::getRegions()], [filter.qseaSet()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   # keep short windows as a toy example
#'   filterWindows(exampleTumourNormal, width < 500)
#' }
#' }
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

#' Arrange (reorder) samples in a qseaSet via dplyr syntax
#'
#' Extends [dplyr::arrange()] to reorder **samples** of a `qseaSet` based on columns
#' in its `sampleTable`.
#'
#' @method arrange qseaSet
#' @importFrom dplyr arrange
#' @param .data A `qseaSet`.
#' @param ...  Variables to arrange by (passed to [dplyr::arrange()]).
#' @param .by_group Ignored; grouping is not used for `qseaSet`.
#'
#' @return The input `qseaSet` object with sample order updated.
#'
#' @seealso [qsea::getSampleTable()], [subsetQset()]
#'
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   exampleTumourNormal %>% arrange(age)
#' }
#' }
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

#' Column names of a qseaSet sample table
#'
#' S4 method for BiocGenerics::colnames() that returns the column names of a
#' qseaSet's sample metadata table (qsea::getSampleTable(x)).
#'
#' @param x A `qseaSet`.
#' @return Character vector of column names.
#' @aliases colnames,qseaSet-method
#' @importFrom BiocGenerics colnames
#' @seealso qsea::getSampleTable, qsea::getSampleNames
#' @examples
#' \donttest{
#' if (system.file("data","exampleTumourNormal.rda",package="mesa") != "") {
#'   data(exampleTumourNormal, package="mesa")
#'   colnames(exampleTumourNormal)
#' }
#' }
#' @export
setMethod("colnames", "qseaSet", function(x) {
  colnames(qsea::getSampleTable(x))
})

