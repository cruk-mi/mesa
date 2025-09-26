#' Pipe operator
#'
#' Re-export of the magrittr pipe. See
#' \code{\link[magrittr:pipe]{\%>\%}} for full details.
#' It forwards the left-hand side (LHS) as the first argument to the
#' right-hand side (RHS) call.
#'
#' @name %>%
#' @rdname pipe
#' @usage lhs \%>\% rhs
#' @return The result of applying \code{rhs} to \code{lhs}.
#' @keywords internal
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Basic usage
#' 1:5 %>% sum()
#'
#' # Chain multiple transformations
#' 1:10 %>%
#'   mean() %>%
#'   round(1)
NULL
