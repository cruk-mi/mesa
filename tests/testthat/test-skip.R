test_that("skip_long_checks() skips only when the option is TRUE", {

    old <- options(skip_long_checks = NULL)
    on.exit(options(old), add = TRUE)

    expect_true(skip_long_checks())

    options(skip_long_checks = FALSE)
    expect_true(skip_long_checks())

    options(skip_long_checks = TRUE)
    expect_condition(skip_long_checks(), class = "skip")
})
