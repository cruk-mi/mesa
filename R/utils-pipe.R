#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' 
#' @examples
#' # Basic usage: pipe a vector into sum()
#' 1:5 %>% sum()
#'
#' # Chain multiple transformations
#' 1:10 %>%
#'   mean() %>%
#'   round(1)
NULL
