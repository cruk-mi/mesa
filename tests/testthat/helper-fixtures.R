# Cache qsea::getExampleQseaSet() results across test files within a
# single R session, so repeated builds at the same (repl, depth) cost
# only one simulation. testthat sources helper-*.R once per test run,
# before any test_that() block.

.mesa_test_fixtures <- new.env(parent = emptyenv())

cachedExampleQset <- function(repl = 8, expSamplingDepth = 1e5) {
    key <- paste0("repl", repl, "_d", format(expSamplingDepth, scientific = FALSE))
    if (is.null(.mesa_test_fixtures[[key]])) {
        .mesa_test_fixtures[[key]] <- qsea::getExampleQseaSet(
            repl = repl, expSamplingDepth = expSamplingDepth
        )
    }
    .mesa_test_fixtures[[key]]
}
