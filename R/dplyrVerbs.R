#' Filter samples in a qseaSet
#'
#' Extends [dplyr::filter()] to subset **samples** in a `qseaSet` based on
#' predicates applied to its `sampleTable`.
#'
#' @details
#' Filtering is performed on the result of `qsea::getSampleTable(.data)`.
#' Matching rows determine the `sample_name`s retained in the `qseaSet`.
#' Grouped filtering (`group_by`) is not supported for `qseaSet`.
#'
#' @method filter qseaSet
#' @importFrom dplyr filter
#'
#' @param .data `qseaSet`  
#'   A `qseaSet` object whose samples are to be filtered.
#'
#' @param ...  
#'   Predicate expressions passed to [dplyr::filter()], evaluated on the
#'   `sampleTable`. For example, `age > 65`, `tissue == "Colon"`.
#'   **Default:** none (no filtering).
#'
#' @param .preserve `logical(1)`  
#'   Placeholder argument for compatibility with `dplyr::filter()`. Ignored
#'   here since grouping is not implemented.  
#'   **Default:** `FALSE`.
#'
#' @return A `qseaSet` object containing only the samples whose
#'   `sampleTable` rows satisfied the filter expressions. If no samples
#'   remain, an empty `qseaSet` is returned (with a message but not an error).
#'
#' @seealso
#' [qsea::getSampleTable()] to inspect metadata,  
#' [subsetQset()] for explicit sample subsetting,  
#' [selectQset()] for column selection in the `sampleTable`.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Keep only older patients
#' exampleTumourNormal %>% filter(age >= 70)
#'
#' # Restrict to colon tissue
#' exampleTumourNormal %>% filter(tissue == "Colon")
#'
#' # Select only tumour samples
#' exampleTumourNormal %>% filter(tumour != "Normal")
#'
#' # Combine predicates: lung tumours only
#' exampleTumourNormal %>%
#'   filter(stringr::str_detect(sample_name, "Lung"), tumour != "Normal")
#'
#' # If no rows match, returns empty qseaSet (with a message)
#' exampleTumourNormal %>% filter(tissue == "Liver")
#'
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
#' Extends [dplyr::mutate()] to add or modify columns in a `qseaSet`'s
#' `sampleTable`.
#'
#' @details
#' The order of samples is preserved.  
#' The `sample_name` column **cannot** be changed with this function; use
#' [renameQsetNames()] or [renameSamples()] if sample names must be updated.
#'
#' @method mutate qseaSet
#' @importFrom dplyr mutate
#'
#' @param .data `qseaSet`  
#'   A `qseaSet` object whose `sampleTable` will be mutated.
#'
#' @param ...   
#'   Column transformations passed to [dplyr::mutate()]. For example,
#'   `over70 = age > 70`, `group = paste(tissue, tumour, sep = "_")`.
#'   **Default:** none.
#'
#' @return A `qseaSet` object with its `sampleTable` updated:  
#' * new columns are added,  
#' * existing columns are modified,  
#' * `sample_name` remains unchanged.
#'
#' @seealso
#' [qsea::getSampleTable()] to inspect metadata,  
#' [renameSamples()] and [renameQsetNames()] for renaming,  
#' [filter.qseaSet()] and [selectQset()] for other dplyr verbs.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Add a new logical column based on age
#' exampleTumourNormal %>%
#'   mutate(over70 = age > 70) %>%
#'   qsea::getSampleTable()
#'
#' # Recode gender (better practice: use dplyr::case_when or stringr::str_replace)
#' exampleTumourNormal %>%
#'   mutate(gender = ifelse(gender == "F", "Female", "Male")) %>%
#'   qsea::getSampleTable()
#'
#' # Change group by removing digits from sample_name
#' exampleTumourNormal %>%
#'   mutate(group = stringr::str_remove(sample_name, "[0-9]")) %>%
#'   qsea::getSampleTable()
#' 
#' @export
mutate.qseaSet <- function(.data, ...){

   newTable <- .data %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column(".rownameCol") %>%
    dplyr::mutate(...) %>%
    tibble::column_to_rownames(".rownameCol")
  
  if (!(identical(.data@sampleTable$sample_name, newTable$sample_name))) {
    stop(glue::glue("sample_name cannot be changed with dplyr::mutate(). Use mesa::renameQsetNames() or mesa::renameSamples() instead."))
  }

  .data@sampleTable <- newTable

  if (!(identical(.data@sampleTable$sample_name, newTable$sample_name))) {
    stop(glue::glue("sample_name cannot be changed with dplyr::mutate(). Use mesa::renameQsetNames() or mesa::renameSamples() instead."))
  }

  .data@sampleTable <- newTable

  return(.data)}


#' Left join data onto a qseaSet sample table
#'
#' Extends [dplyr::left_join()] to merge a data frame into a `qseaSet`'s
#' `sampleTable`. Useful for adding annotations or metadata to samples.
#'
#' @method left_join qseaSet
#' @importFrom dplyr left_join
#'
#' @param x `qseaSet`.  
#'   A qseaSet object whose `sampleTable` will receive additional columns.
#'
#' @param y `data.frame`.  
#'   A data frame or tibble to join onto the `sampleTable`.
#'
#' @param by `character()` or `named character()`.  
#'   Join keys, as in [dplyr::left_join()].  
#'   For qseaSets, the safest choice is `by = "sample_name"`.  
#'   **Default:** `NULL` (all common variables between `x@sampleTable` and `y`).
#'
#' @param copy `logical(1)`.  
#'   For remote backends only. If `TRUE`, copies `y` into the same data source
#'   as `x`. Ignored for in-memory data frames.  
#'   **Default:** `FALSE`.
#'
#' @param suffix `character(2)`.  
#'   Suffixes appended to overlapping non-join column names from `x` and `y`.  
#'   **Default:** `c(".x", ".y")`.
#'
#' @param keep `logical(1)` or `NULL`.  
#'   If `TRUE`, keeps join keys from both `x` and `y`.  
#'   If `FALSE`, only keys from `x` are kept.  
#'   **Default:** `NULL` (deferred to [dplyr::left_join()] behaviour).
#'
#' @param ... Additional arguments passed to [dplyr::left_join()],
#'   e.g. `relationship = "one-to-one"` (dplyr ≥ 1.1).
#'   **Default:** none.
#'
#' @return A `qseaSet` object with its `sampleTable` updated:  
#' * existing rows preserved,  
#' * additional columns from `y` merged by sample name or other keys,  
#' * duplicate column names disambiguated with `suffix`.  
#'
#' @seealso
#' [qsea::getSampleTable()] to view metadata,  
#' [mutate.qseaSet()], [filter.qseaSet()] for other dplyr verbs,  
#' [dplyr::left_join()] for join semantics.
#'
#' @examples
#'  data(exampleTumourNormal, package="mesa")
#'   
#'  # Join 'patient' (y) to 'sample_name' (x)
#'  newData <- tibble::tibble(patient = c("Colon1","Colon2","Lung1","Lung3"), new = 1:4)
#'  exampleTumourNormal %>%
#'   left_join(newData, by = c("sample_name" = "patient")) %>%
#'   qsea::getSampleTable()
#' 
#' @export
left_join.qseaSet <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x",".y"), keep = NULL, ...){

  x@sampleTable <- x %>%
    qsea::getSampleTable() %>%
    tibble::rownames_to_column("rownameCol") %>%
    dplyr::left_join(y, by = by, copy = copy, suffix = suffix, keep = keep,...) %>%
    tibble::column_to_rownames("rownameCol")

  return(x)}


#' Select or rename columns in a qseaSet sample table
#'
#' Extends [dplyr::select()] to keep or rename columns in a `qseaSet`'s
#' `sampleTable`. This allows column subsetting, renaming, and use of
#' tidyselect helpers while preserving key identifiers.
#'
#' @details
#' This is a wrapper around [selectQset()].  
#' Essential columns `sample_name` and `group` are always preserved (and appear
#' first if present).
#'
#' @method select qseaSet
#' @importFrom dplyr select
#'
#' @param .data `qseaSet`.  
#'   A qseaSet object whose `sampleTable` columns will be selected or renamed.
#'
#' @param ...
#'   Additional arguments passed to [dplyr::select()], supports tidyselect helpers
#'   (e.g., [dplyr::matches()], [dplyr::starts_with()]).  
#'   Negative selection (`-col`) is supported, but `sample_name` and `group`
#'   cannot be removed.
#'   **Default:** none.
#'
#' @return A `qseaSet` with its `sampleTable` updated:  
#' * contains only selected/renamed columns,  
#' * always includes `sample_name` and `group` (if available).  
#'
#' @seealso
#' [selectQset()] for the underlying implementation,  
#' [qsea::getSampleTable()] to view results,  
#' [mutate.qseaSet()], [filter.qseaSet()] for related dplyr verbs.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Keep only specific columns (sample_name and group always preserved)
#' exampleTumourNormal %>%
#'   dplyr::select(type, tumour) %>%
#'   qsea::getSampleTable()
#'
#' # Select and rename a column
#' exampleTumourNormal %>%
#'   dplyr::select(age, tumour_new = tumour) %>%
#'   qsea::getSampleTable()
#'
#' # Use tidyselect helpers
#' exampleTumourNormal %>%
#'   dplyr::select(matches("s")) %>%
#'   qsea::getSampleTable()
#'
#' # Negative selection (except essential columns)
#' exampleTumourNormal %>%
#'   dplyr::select(-matches("s")) %>%
#'   qsea::getSampleTable()
#'
#' @export
select.qseaSet <- function(.data, ...){selectQset(.data, ...)}


#' Select or rename columns of a qseaSet sample table
#'
#' Keeps or renames columns from the `sampleTable` of a `qseaSet`.  
#' Essential columns `sample_name` and `group` are always preserved (and placed
#' at the front of the table if present).
#'
#' @param qseaSet `qseaSet`.  
#'   A qseaSet object whose `sampleTable` columns will be selected or renamed.
#'
#' @param ...  
#'   Additional arguments passed to [dplyr::select()]. Supports tidyselect helpers
#'   (e.g., [dplyr::matches()], [dplyr::starts_with()]) and renaming
#'   (`newName = oldName`). Negative selection is supported, but `sample_name`
#'   and `group` cannot be removed.
#'   **Default:** none.
#'
#' @return A `qseaSet` with its `sampleTable` updated:  
#' * contains only the selected/renamed columns,  
#' * always includes `sample_name` and `group` (if available).  
#'
#' @seealso
#' [select.qseaSet()] for the dplyr S3 method,  
#' [qsea::getSampleTable()] to inspect the updated table.
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Keep only sample_name and group
#' exampleTumourNormal %>%
#'   selectQset(sample_name, group) %>%
#'   qsea::getSampleTable()
#'
#' # Select and rename a column
#' exampleTumourNormal %>%
#'   selectQset(age, tumour_type = tumour) %>%
#'   qsea::getSampleTable()
#'
#' # Use tidyselect helpers
#' exampleTumourNormal %>%
#'   selectQset(matches("s")) %>%
#'   qsea::getSampleTable()
#'
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
#' Extends [dplyr::pull()] to extract a column from the `sampleTable`
#' of a `qseaSet`.
#'
#' @method pull qseaSet
#' @importFrom dplyr pull
#'
#' @param .data `qseaSet`.  
#'   A qseaSet object whose `sampleTable` column should be extracted.
#'
#' @param var Column to extract (name or position).  
#'   **Default:** `-1` (last column).
#'
#' @param name Optional column to use for names in the returned vector.  
#'   **Default:** `NULL` (no names).
#'
#' @param ... Additional arguments passed to [dplyr::pull()].
#'   **Default:** none.
#'
#' @return A vector with length equal to the number of samples in the `qseaSet`.  
#' If `name` is provided, the vector is named accordingly.
#'
#' @seealso
#' [qsea::getSampleTable()],  
#' [dplyr::pull()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Pull the 'group' column as a vector
#' exampleTumourNormal %>% pull(group)
#'
#' # Pull by column position (last column by default)
#' exampleTumourNormal %>% pull()
#'
#' # Pull a column and name it by sample_name
#' exampleTumourNormal %>% pull(group, name = sample_name)
#'
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
#' Extends [base::sort()] to reorder the **samples** of a `qseaSet`
#' using [gtools::mixedsort()] so that numeric parts sort naturally
#' (e.g., `"2"` comes before `"10"`).
#'
#' @method sort qseaSet
#'
#' @param x `qseaSet`.  
#'   A qseaSet object to reorder.
#'
#' @param decreasing `logical(1)`.  
#'   Whether to reverse the order.  
#'   **Default:** `FALSE`.
#'
#' @param ... Additional arguments passed on to [gtools::mixedsort()].  
#'   **Default:** none.
#'
#' @return A `qseaSet` object with reordered samples. The content of all slots
#' (e.g., `count_matrix`, `sampleTable`) is updated consistently.
#'
#' @seealso
#' [qsea::getSampleNames()],  
#' [subsetQset()],  
#' [gtools::mixedsort()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Sort samples in natural (mixed) order
#' exampleTumourNormal %>% sort()
#'
#' # Sort in reverse order
#' exampleTumourNormal %>% sort(decreasing = TRUE)
#'
#' @export
sort.qseaSet <- function(x, decreasing = FALSE, ...){

  x %>%
    subsetQset(samplesToKeep = (x %>% qsea::getSampleNames() %>% gtools::mixedsort(decreasing = decreasing))) %>%
    return()

  }


#' Filter regions (windows) inside a qseaSet
#'
#' Filters the **regions** of a `qseaSet` using [dplyr::filter()] predicates
#' applied to its regions (as a data frame), then subsets the object accordingly.
#'
#' @param qseaSet `qseaSet`.  
#'   The qseaSet object whose regions (windows) will be filtered.
#'
#' @param ... Expressions.  
#'   Predicates passed to [dplyr::filter()] and evaluated on the regions data frame.  
#'   **Default:** none (no filtering).
#'
#' @return A `qseaSet` object with only regions matching the filter conditions.  
#' Both the `count_matrix` and `@regions` slot are updated consistently.
#'
#' @seealso
#' [qsea::getRegions()],  
#' [filter.qseaSet()] (for filtering samples instead of regions),  
#' [dplyr::filter()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Keep regions starting before position 25,020,000
#' exampleTumourNormal %>% filterWindows(start <= 25020000)
#'
#' # Keep regions with CpG density > 10
#' exampleTumourNormal %>% filterWindows(CpG_density > 10)
#'
#' # Combine multiple conditions
#' exampleTumourNormal %>%
#'   filterWindows(seqnames == "chr7", CpG_density > 5)
#' 
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
#'
#' @param .data `qseaSet`.  
#'   The qseaSet whose samples will be reordered.
#'
#' @param ... Expressions.  
#'   Variables or helper functions passed to [dplyr::arrange()] to define the
#'   ordering.  
#'   **Default:** none (no reordering).
#'
#' @param .by_group `logical(1)`.  
#'   Grouped arrangement is not implemented for `qseaSet`; this argument is ignored.  
#'   **Default:** `FALSE`.
#'
#' @return A `qseaSet` object with samples reordered.  
#' Both the `count_matrix` and the `sampleTable` are updated consistently.
#'
#' @seealso
#' [qsea::getSampleTable()],  
#' [subsetQset()],  
#' [sort.qseaSet()] (for natural sorting of sample names).
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Reorder samples by age
#' exampleTumourNormal %>% arrange(age)
#'
#' # Reorder by tissue
#' exampleTumourNormal %>% arrange(tissue)
#'
#' # Reorder by tissue (descending), then tumour type
#' exampleTumourNormal %>% arrange(desc(tissue), tumour)
#'
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
#' S4 method for \code{\link[BiocGenerics]{colnames}} that returns the column names of a
#' `qseaSet`'s sample metadata table (i.e., [qsea::getSampleTable()]).
#'
#' @param x `qseaSet`.  
#'   The qseaSet object whose `sampleTable` column names are to be retrieved.
#'
#' @return `character()`.  
#'   A character vector containing the column names of the sample table.
#'
#' @aliases colnames,qseaSet-method
#' @importFrom BiocGenerics colnames
#'
#' @seealso
#' [qsea::getSampleTable()],  
#' [qsea::getSampleNames()]
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Retrieve the column names of the sample metadata
#' colnames(exampleTumourNormal)
#'
#' @exportMethod colnames
setMethod("colnames", "qseaSet", function(x) {
  colnames(qsea::getSampleTable(x))
})

