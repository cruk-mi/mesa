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

test_that("mesaPCA and mesaUMAP accessors return their slots", {
    mp <- validPCA()
    expect_identical(getPrcomp(mp), makePrcomp())
    expect_identical(getCoordinates(mp), as.data.frame(makePrcomp()$x))
    expect_identical(getWindowNames(mp), paste0("W", 1:4))

    mu <- validUMAP()
    expect_identical(getCoordinates(mu), makePoints())
    expect_identical(getWindowNames(mu), paste0("W", 1:4))
})

test_that("mesaDimRed accessors return their slots", {
    res <- list(pca1 = validPCA(), pca2 = validPCA())
    dt <- data.frame(S1 = 0, S2 = 0, S3 = 0, S4 = 0, S5 = 0)
    params <- list(method = "PCA", normMethod = "nrpm")
    md <- mesaDimRed(res = res, sampleTable = makeSampleTable(),
        params = params, dataTable = dt)

    expect_identical(getResults(md), res)
    expect_identical(getParameters(md), params)
    expect_identical(getDimRedData(md), dt)
    expect_identical(getCoordinates(md),
        list(pca1 = getCoordinates(res$pca1),
            pca2 = getCoordinates(res$pca2)))
    expect_identical(getWindowNames(md),
        list(pca1 = paste0("W", 1:4), pca2 = paste0("W", 1:4)))

    # An empty container gives empty results.
    empty <- mesaDimRed(res = list(), sampleTable = makeSampleTable(),
        params = list())
    expect_identical(getResults(empty), list())
    expect_identical(getDimRedData(empty), data.frame())
    expect_length(getCoordinates(empty), 0)
    expect_length(getWindowNames(empty), 0)
})

test_that("as.data.frame() gives one long table of every result", {
    st <- makeSampleTable()
    md <- mesaDimRed(res = list(pca1 = validPCA(), pca2 = validPCA()),
        sampleTable = st, params = list(method = "PCA"))
    df <- as.data.frame(md)

    expect_s3_class(df, "data.frame")
    expect_identical(nrow(df), 10L)
    expect_identical(df$resName, rep(c("pca1", "pca2"), each = 5))
    expect_identical(df$sample_name, rep(paste0("S", 1:5), 2))
    pcs <- colnames(getCoordinates(validPCA()))
    expect_identical(colnames(df), c("resName", "sample_name", pcs, "group"))
    expect_equal(as.matrix(df[1:5, pcs]), makePrcomp()$x,
        ignore_attr = TRUE)
    expect_identical(df$group, rep(st$group, 2))

    umap <- mesaDimRed(res = list(umap1 = validUMAP()), sampleTable = st,
        params = list(method = "UMAP"))
    expect_identical(colnames(as.data.frame(umap)),
        c("resName", "sample_name", "UMAP1", "UMAP2", "group"))
})

test_that("as.data.frame() handles group means and empty containers", {
    groups <- mesaUMAP(makePoints(c("A", "B")), "W1")
    grouped <- mesaDimRed(res = list(umap1 = groups),
        sampleTable = makeSampleTable(),
        params = list(method = "UMAP", useGroupMeans = TRUE))
    df <- as.data.frame(grouped)
    expect_identical(colnames(df), c("resName", "group", "UMAP1", "UMAP2"))
    expect_identical(df$group, c("A", "B"))

    empty <- mesaDimRed(res = list(), sampleTable = makeSampleTable(),
        params = list())
    expect_identical(as.data.frame(empty),
        data.frame(resName = character(), sample_name = character()))
})

test_that("getWindowNames() still labels qseaSet, GRanges and data.frame", {
    gr <- GenomicRanges::GRanges(c("chr1", "chr2"),
        IRanges::IRanges(c(10, 20), c(15, 30)))
    expect_identical(getWindowNames(gr), c("chr1:10-15", "chr2:20-30"))
    df <- data.frame(seqnames = c("chr1", "chr2"), start = c(10, 20),
        end = c(15, 30))
    expect_identical(getWindowNames(df), c("chr1:10-15", "chr2:20-30"))
    expect_length(getWindowNames(exampleTumourNormal),
        length(qsea::getRegions(exampleTumourNormal)))
})

test_that("getPCA() accessors agree with the PCA it ran", {
    pca <- getPCA(exampleTumourNormal, topVarNum = c(10, 100),
        verbose = FALSE)
    expect_named(getResults(pca), names(getCoordinates(pca)))
    expect_identical(getParameters(pca)$method, "PCA")
    expect_identical(lengths(getWindowNames(pca), use.names = FALSE),
        c(10L, 100L))
    expect_identical(rownames(getCoordinates(pca)[[1]]),
        getSampleNames(pca))

    df <- as.data.frame(pca)
    expect_identical(nrow(df), 2L * length(getSampleNames(pca)))
    expect_true(all(colnames(getSampleTable(pca)) %in% colnames(df)))
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
