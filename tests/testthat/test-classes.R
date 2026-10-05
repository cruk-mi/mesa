# Small hand-built results, so every validity rule can be broken one at a time.
makePrcomp <- function(samples = paste0("S", 1:5), windows = paste0("W", 1:4)) {
    set.seed(1)
    x <- matrix(rnorm(length(samples) * length(windows)),
        nrow = length(samples), dimnames = list(samples, windows))
    stats::prcomp(x)
}

makePoints <- function(samples = paste0("S", 1:5)) {
    data.frame(UMAP1 = seq_along(samples), UMAP2 = rev(seq_along(samples)),
        row.names = samples)
}

makeSampleTable <- function(samples = paste0("S", 1:5)) {
    data.frame(sample_name = samples, group = rep(c("A", "B"), c(3, 2)),
        row.names = samples)
}

validPCA <- function() {
    mesaPCA(prcomp = makePrcomp(), windows = paste0("W", 1:4))
}

validUMAP <- function() {
    mesaUMAP(points = makePoints(), windows = paste0("W", 1:4))
}

test_that("mesaPCA: a valid object constructs, validates and shows", {
    mp <- validPCA()
    expect_s4_class(mp, "mesaPCA")
    expect_true(validObject(mp))
    expect_output(show(mp),
        "PCA result for 5 samples calculated over 4 windows")
})

test_that("mesaPCA: each validity rule rejects an inconsistent object", {
    pc <- makePrcomp()
    expect_error(mesaPCA(pc, character()), "`windows` must not be empty")
    expect_error(mesaPCA(pc, c("W1", "W2", "W3", NA)),
        "`windows` must not contain NA")
    expect_error(mesaPCA(pc, c("W1", "W2", "W3", "W3")),
        "`windows` must not contain duplicates")

    noX <- pc
    noX$x <- NULL
    expect_error(mesaPCA(noX, paste0("W", 1:4)),
        "`prcomp\\$x` must be a matrix with sample IDs as row names")
    noNames <- pc
    rownames(noNames$x) <- NULL
    expect_error(mesaPCA(noNames, paste0("W", 1:4)),
        "`prcomp\\$x` must be a matrix with sample IDs as row names")

    expect_error(mesaPCA(pc, paste0("W", 1:3)),
        "`windows` must have one entry per row of `prcomp\\$rotation`")
})

test_that("mesaUMAP: a valid object constructs, validates and shows", {
    mu <- validUMAP()
    expect_s4_class(mu, "mesaUMAP")
    expect_true(validObject(mu))
    expect_output(show(mu),
        "UMAP result for 5 samples calculated over 4 windows")
})

test_that("mesaUMAP: each validity rule rejects an inconsistent object", {
    pts <- makePoints()
    expect_error(mesaUMAP(pts, character()), "`windows` must not be empty")
    expect_error(mesaUMAP(pts, c("W1", NA)), "`windows` must not contain NA")
    expect_error(mesaUMAP(pts, c("W1", "W1")),
        "`windows` must not contain duplicates")

    expect_error(mesaUMAP(pts[0, ], "W1"),
        "`points` must have at least one row")
    expect_error(mesaUMAP(data.frame(UMAP1 = 1:3, UMAP2 = 3:1), "W1"),
        "`points` must have sample IDs as row names")
    expect_error(
        mesaUMAP(data.frame(UMAP1 = c("a", "b"), row.names = c("S1", "S2")),
            "W1"),
        "`points` columns must all be numeric"
    )
})

test_that("mesaDimRed: a valid object constructs, validates and shows", {
    md <- mesaDimRed(res = list(pca1 = validPCA()),
        sampleTable = makeSampleTable(), params = list(method = "PCA"))
    expect_s4_class(md, "mesaDimRed")
    expect_true(validObject(md))
    expect_output(show(md),
        "1 dimensionality reduction objects for 5 samples")
    expect_identical(getSampleTable(md), makeSampleTable())
    expect_identical(getSampleNames(md), paste0("S", 1:5))

    # An empty container is valid.
    empty <- mesaDimRed(res = list(), sampleTable = makeSampleTable(),
        params = list())
    expect_true(validObject(empty))
    expect_output(show(empty),
        "0 dimensionality reduction objects for 0 samples")
    expect_identical(getSampleNames(empty), character())
})

test_that("mesaDimRed: each validity rule rejects an inconsistent object", {
    st <- makeSampleTable()
    build <- function(res, params = list(), sampleTable = st,
        dataTable = data.frame()) {
        mesaDimRed(res = res, sampleTable = sampleTable, params = params,
            dataTable = dataTable)
    }

    expect_error(build(list(pca1 = makePrcomp())),
        "every element of `res` must be a mesaPCA or mesaUMAP object")
    expect_error(build(list(pca1 = validPCA(), umap1 = validUMAP())),
        "`res` must not mix mesaPCA and mesaUMAP objects")

    expect_error(build(list(validPCA())),
        "`res` must be a named list with unique names")
    expect_error(build(list(pca1 = validPCA(), validPCA())),
        "`res` must be a named list with unique names")
    expect_error(build(list(pca1 = validPCA(), pca1 = validPCA())),
        "`res` must be a named list with unique names")

    expect_error(build(list(pca1 = validPCA()), params = list(method = "UMAP")),
        "`params\\$method` is \"UMAP\" but `res` holds PCA results")

    reordered <- mesaPCA(makePrcomp(samples = paste0("S", 5:1)),
        paste0("W", 1:4))
    expect_error(build(list(pca1 = validPCA(), pca2 = reordered)),
        "every element of `res` must cover the same samples, in the same order")

    dupSamples <- mesaPCA(makePrcomp(samples = c("S1", "S1", "S2", "S3")),
        paste0("W", 1:4))
    expect_error(build(list(pca1 = dupSamples)),
        "sample IDs in `res` must be unique")

    unknown <- mesaPCA(makePrcomp(samples = c(paste0("S", 1:4), "S9")),
        paste0("W", 1:4))
    expect_error(build(list(pca1 = unknown)),
        "sample IDs in `res` are missing from `sampleTable`: S9")

    # With group means the IDs are groups, checked against sampleTable$group.
    groups <- mesaUMAP(makePoints(c("A", "B")), "W1")
    grouped <- build(list(umap1 = groups),
        params = list(method = "UMAP", useGroupMeans = TRUE))
    expect_true(validObject(grouped))
    expect_identical(getSampleNames(grouped), c("A", "B"))
    expect_error(build(list(umap1 = groups), params = list(method = "UMAP")),
        "sample IDs in `res` are missing from `sampleTable`: A, B")

    dt <- data.frame(seqnames = "1", start = 1, end = 2,
        S1 = 0, S2 = 0, S3 = 0, S4 = 0)
    expect_error(build(list(pca1 = validPCA()), dataTable = dt),
        "`dataTable` has no column for samples: S5")
    expect_true(validObject(build(list(pca1 = validPCA()),
        dataTable = cbind(dt, S5 = 0))))
})

test_that("getPCA() and getUMAP() output passes validObject()", {
    qs <- cachedExampleQset()
    expect_true(validObject(getPCA(qs, normMethod = "nrpm", verbose = FALSE)))
    expect_identical(
        getSampleNames(getPCA(exampleTumourNormal, verbose = FALSE)),
        qsea::getSampleNames(exampleTumourNormal)
    )

    expect_true(validObject(getPCA(exampleTumourNormal,
        topVarNum = c(10, 100), returnDataTable = TRUE, verbose = FALSE)))

    expect_true(validObject(getPCA(exampleTumourNormal,
        dataTable = getBetaTable(exampleTumourNormal), verbose = FALSE)))

    grouped <- exampleTumourNormal %>%
        mutate(group = stringr::str_remove(sample_name, "_[NT]$"))
    expect_true(validObject(getPCA(grouped, useGroupMeans = TRUE,
        returnDataTable = TRUE, verbose = FALSE)))

    set.seed(1)
    expect_true(validObject(getUMAP(exampleTumourNormal, n_neighbors = 5,
        returnDataTable = TRUE)))
})
