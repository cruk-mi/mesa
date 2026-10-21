# Skip a slow test when options(skip_long_checks = TRUE) is set.
skip_long_checks <- function() {
    if (!isTRUE(getOption("skip_long_checks"))) {
        return(invisible(TRUE))
    }

    testthat::skip(
        "Slow checks skipped when options(skip_long_checks = TRUE) has been set"
    )
}
