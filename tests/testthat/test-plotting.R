test_that("ensemblReachable() is TRUE when the BioMart service answers", {

    local_mocked_bindings(martServiceStatus = function(url, timeout) 200L)

    expect_true(ensemblReachable("jul2022.archive.ensembl.org"))
})

test_that("ensemblReachable() is FALSE on an HTTP error or a timeout", {

    local_mocked_bindings(martServiceStatus = function(url, timeout) 500L)
    expect_false(ensemblReachable("jul2022.archive.ensembl.org"))

    # curlGetHeaders() signals a timeout as an error (#183)
    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        stop("libcurl error code 28:\n\tOperation timed out")
    })
    expect_false(ensemblReachable("jul2022.archive.ensembl.org"))
})

test_that("ensemblReachable() probes the BioMart service, not the site root", {

    seen <- NULL
    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        seen <<- list(url = url, timeout = timeout)
        200L
    })

    # the jul2023 root answers while its BioMart service does not (#183)
    service <- "/biomart/martservice?type=registry"
    ensemblReachable("jul2023.archive.ensembl.org", timeout = 2)
    expect_equal(seen$url,
        paste0("https://jul2023.archive.ensembl.org", service))
    expect_equal(seen$timeout, 2)

    ensemblReachable("https://jul2022.archive.ensembl.org")
    expect_equal(seen$url,
        paste0("https://jul2022.archive.ensembl.org", service))
    expect_equal(seen$timeout, 5)
})

test_that("ensemblReachable() returns FALSE at once when the probe errors", {

    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        stop("Connection timed out")
    })

    elapsed <- system.time(reachable <- ensemblReachable("10.255.255.1"))
    expect_false(reachable)
    expect_lt(elapsed[["elapsed"]], 1)
})

test_that("ensemblReachable() rounds the timeout up to whole seconds", {

    # curlGetHeaders() reads 0 as its 60 s default, so 0.5 must not become 0
    sent <- NULL
    local_mocked_bindings(martServiceStatus = function(url, timeout) {
        sent <<- timeout
        200L
    })

    ensemblReachable("jul2022.archive.ensembl.org", timeout = 0.5)
    expect_equal(sent, 1)

    ensemblReachable("jul2022.archive.ensembl.org", timeout = 2.2)
    expect_equal(sent, 3)
})

test_that("ensemblReachable() rejects a timeout that is not positive", {

    local_mocked_bindings(martServiceStatus = function(url, timeout) 200L)

    for (bad in list(0, -1, NA_real_, "5", c(1, 2))) {
        expect_error(
            ensemblReachable("jul2022.archive.ensembl.org", timeout = bad),
            "positive number"
        )
    }
})

test_that("plotGeneHeatmap() fails fast when the default archive hangs", {

    # purrr::possibly() swallows errors from useMart(), so record the call
    called <- FALSE
    local_mocked_bindings(ensemblReachable = function(host, ...) FALSE)
    local_mocked_bindings(
        useMart = function(...) {
            called <<- TRUE
            stop("useMart() must not be called")
        },
        .package = "biomaRt"
    )

    expect_error(
        plotGeneHeatmap(exampleTumourNormal, gene = "HOXA9"),
        "did not answer within 5 s"
    )
    expect_false(called)
})

test_that("plotGeneHeatmap() waits for mesa_ensemblTimeout seconds", {

    oldOptions <- options(mesa_ensemblTimeout = 12)
    on.exit(options(oldOptions), add = TRUE)

    waited <- NULL
    local_mocked_bindings(ensemblReachable = function(host, timeout) {
        waited <<- timeout
        FALSE
    })

    expect_error(
        plotGeneHeatmap(exampleTumourNormal, gene = "HOXA9"),
        "did not answer within 12 s"
    )
    expect_equal(waited, 12)
})
