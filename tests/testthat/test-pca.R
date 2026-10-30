test_that("PCAs", {

    library(rlang)

    # Test that getPCA() returns a valid result for various input arguments
    expect_no_error(obj1 <- exampleTumourNormal %>% getPCA())
    expect_no_error(obj2 <- exampleTumourNormal %>% getPCA(returnDataTable = TRUE))
    expect_no_error(obj3 <- exampleTumourNormal %>% getPCA(topVarNum = 10))
    expect_no_error(obj4 <- exampleTumourNormal %>% getPCA(topVarSamples = "_T", topVarNum = 10))
    expect_no_error(obj5 <- exampleTumourNormal %>% getPCA(minDensity = 10))
    expect_no_error(obj6 <- exampleTumourNormal %>% getPCA(topVarNum = c(10,100,200)))
    expect_no_error(obj7 <- exampleTumourNormal %>% getPCA(dataTable = getBetaTable(exampleTumourNormal)))

    expect_error(exampleTumourNormal %>% filter(str_detect(sample_name, "Colon1_T")) %>% getPCA())
    expect_error(exampleTumourNormal %>% filter(str_detect(sample_name, "Colon1")) %>% getPCA())
    expect_error(exampleTumourNormal %>% filterWindows(seqnames == 1) %>% getPCA()) #no windows left

    # Test that the x matrices in the res components are not empty
    expect_false(is_empty(getCoordinates(obj1)$pca1))
    expect_false(is_empty(getCoordinates(obj2)$pca1))
    expect_false(is_empty(getCoordinates(obj3)$pca1))
    expect_false(is_empty(getCoordinates(obj4)$pca1))
    expect_false(is_empty(getCoordinates(obj5)$pca1))
    expect_false(is_empty(getCoordinates(obj6)$pca1))
    expect_false(is_empty(getCoordinates(obj7)$pca1))

    # Test that the dimensions of the x matrices are consistent
    expect_equal(dim(getCoordinates(obj1)$pca1), dim(getCoordinates(obj2)$pca1))
    expect_equal(dim(getCoordinates(obj1)$pca1), dim(getCoordinates(obj3)$pca1))
    expect_equal(dim(getCoordinates(obj1)$pca1), dim(getCoordinates(obj4)$pca1))
    expect_equal(dim(getCoordinates(obj1)$pca1), dim(getCoordinates(obj5)$pca1))
    expect_equal(dim(getCoordinates(obj1)$pca1), dim(getCoordinates(obj6)$pca1))
    expect_equal(dim(getCoordinates(obj1)$pca1), dim(getCoordinates(obj7)$pca1))

    # Test that the variance_explained component is a numeric vector
    expect_true(is.numeric(getPrcomp(getResults(obj1)$pca1)$sdev))

    # Test that getPCA() returns an error when called with invalid arguments
    expect_error(exampleTumourNormal %>% getPCA("invalid_argument"))

    # Check res components are of correct length
    expect_equal(length(getResults(obj1)), 1)
    expect_equal(length(getResults(obj6)), 3)

    # Check PCA results are different for different numbers of windows
    expect_false(isTRUE(all.equal(
        getCoordinates(obj6)$pca1,
        getCoordinates(obj6)$pca2)))
    expect_false(isTRUE(all.equal(
        getCoordinates(obj6)$pca1,
        getCoordinates(obj6)$pca3)))

    # Check dimensions of x are correct
    expect_equal(dim(getCoordinates(obj6)$pca1), c(10, 5))
    expect_equal(dim(getCoordinates(obj6)$pca2), c(10, 5))
    expect_equal(dim(getCoordinates(obj6)$pca3), c(10, 5))

    # Check dimensions of rotation are correct
    expect_equal(dim(getPrcomp(getResults(obj6)$pca1)$rotation), c(10, 5))
    expect_equal(dim(getPrcomp(getResults(obj6)$pca2)$rotation), c(100, 5))
    expect_equal(dim(getPrcomp(getResults(obj6)$pca3)$rotation), c(200, 5))

    ##### plotting tests

    # Check no error is returned
    expect_no_error(plots1 <- plotPCA(obj1))
    expect_no_error(plots6 <- plotPCA(obj6))

    # Check output is correct length
    expect_equal(length(plots1), 1)
    expect_equal(length(plots6), 3)

    # Test that plotPCA() returns a ggplot object
    expect_true( "ggplot" %in% class(plots1[[1]][[1]]))
    expect_true( "ggplot" %in% class(plots6[[1]][[1]]))

    # Test other arguments of plotPCA()
    expect_no_error(plotPCA(obj1, colour = "type"))
    expect_no_error(plotPCA(obj1, colour = "age"))
    expect_no_error(plotPCA(obj1, colour = "gender"))
    expect_no_error(plotPCA(obj1 %>% mutate(diver = seq(-4,5)), colour = "diver"))
    # A sampleTable column named like a component must not break the plot
    clash <- plotPCA(obj1 %>% mutate(PC1 = 0))
    expect_no_error(ggplot2::ggplot_build(clash[[1]][[1]]))
    expect_no_error(plotPCA(obj1, shape = "type"))

    expect_no_error(plotPCA(obj1, colour = "type", colourPalette = RColorBrewer::brewer.pal(5,"Oranges")))
    expect_no_error(plotPCA(obj1, colour = "gender", shapePalette = c(2,5), shape = "gender"))

    expect_no_error(plotPCA(obj1 %>% mutate(gender = as.factor(gender)),
                            colour = "gender", shapePalette = c(2,5), shape = "gender"))

    expect_no_error(obj1 %>% mutate(newCol = rnorm(10)) %>% plotPCA(colour = "newCol"))

}

)

test_that("UMAPs", {

    # uwot calls gc() four times per fit, about 0.35 s each with mesa loaded
    local_mocked_bindings(gc = function(...) invisible(NULL), .package = "base")
    set.seed(1)
    expect_no_error(obj1 <- exampleTumourNormal %>% getUMAP(n_neighbors = 5) )
    set.seed(1)
    expect_no_error(obj2 <- exampleTumourNormal %>% getUMAP(returnDataTable = TRUE, n_neighbors = 5) )
    set.seed(1)
    expect_no_error(obj3 <- exampleTumourNormal %>% getUMAP(topVarNum = 10, n_neighbors = 5) )
    set.seed(1)
    expect_no_error(obj4 <- exampleTumourNormal %>% getUMAP(topVarSamples = "_T", topVarNum = 10, n_neighbors = 5) )
    set.seed(1)
    expect_no_error(obj5 <- exampleTumourNormal %>% getUMAP(minDensity = 10, n_neighbors = 5) )
    set.seed(1)
    expect_no_error(obj6 <- exampleTumourNormal %>% getUMAP(topVarNum = c(10,100,200), n_neighbors = 5) )
    set.seed(1)
    expect_no_error(obj7 <- exampleTumourNormal %>% getUMAP(dataTable = getBetaTable(exampleTumourNormal), n_neighbors = 5) )

    expect_error(exampleTumourNormal %>% filter(str_detect(sample_name, "Colon1_T")) %>% getUMAP())
    expect_error(exampleTumourNormal %>% filter(str_detect(sample_name, "Colon1")) %>% getUMAP())
    expect_error(exampleTumourNormal %>% filterWindows(seqnames == 1) %>% getUMAP()) #no windows left

    expect_true(isTRUE(all.equal(
            getCoordinates(obj1)$umap1,
            getCoordinates(obj2)$umap1)))
    expect_false(isTRUE(all.equal(
            getCoordinates(obj3)$umap1,
            getCoordinates(obj4)$umap1)))
    expect_false(isTRUE(all.equal(
            getCoordinates(obj1)$umap1,
            getCoordinates(obj5)$umap1)))
    expect_false(isTRUE(all.equal(
            getCoordinates(obj1)$umap1,
            getCoordinates(obj6)$umap1)))
    expect_true(isTRUE(all.equal(
            getCoordinates(obj1)$umap1,
            getCoordinates(obj7)$umap1)))

    expect_equal(length(getResults(obj1)), 1)

    expect_equal(length(getResults(obj6)), 3)
    expect_false(isTRUE(all.equal(
            getCoordinates(obj6)$umap1,
            getCoordinates(obj6)$umap2)))
    expect_false(isTRUE(all.equal(
            getCoordinates(obj6)$umap1,
            getCoordinates(obj6)$umap3)))

    expect_no_error(plots1 <- plotUMAP(obj1))
    expect_no_error(plots6 <- plotUMAP(obj6))

    expect_no_error(plotUMAP(obj1, colour = "type"))
    expect_no_error(plotUMAP(obj1, colour = "age"))
    expect_no_error(plotUMAP(obj1, colour = "gender"))
    expect_no_error(plotUMAP(obj1 %>% mutate(diver = seq(-4,5)), colour = "diver"))
    expect_no_error(plotUMAP(obj1, shape = "type"))

    expect_no_error(plotUMAP(obj1, colour = "type", colourPalette = RColorBrewer::brewer.pal(5,"Oranges")))
    expect_no_error(plotUMAP(obj1, colour = "gender", shapePalette = c(2,5), shape = "gender"))

    expect_no_error(obj1 %>% mutate(newCol = rnorm(10)) %>% plotUMAP(colour = "newCol"))

    expect_equal(length(plots1), 1)
    expect_equal(length(plots6), 3)

}

)

test_that("show() and plot labels report the samples in each result", {

    # uwot calls gc() four times per fit, about 0.35 s each with mesa loaded
    local_mocked_bindings(gc = function(...) invisible(NULL), .package = "base")

    plotLabels <- function(obj) {
        p <- plotDimRed(obj)[[1]][[1]]
        c(p$labels$title, p$labels$subtitle)
    }

    pcaAll <- getPCA(exampleTumourNormal, verbose = FALSE)
    expect_output(show(pcaAll),
        "1 dimensionality reduction objects for 10 samples")
    expect_output(show(getResults(pcaAll)$pca1),
        "PCA result for 10 samples calculated over 599 windows")
    expect_equal(plotLabels(pcaAll),
        c("PCA for 10 samples using all 599 windows.", "Using beta values."))

    pcaTop <- getPCA(exampleTumourNormal, topVarNum = 10, verbose = FALSE)
    expect_equal(plotLabels(pcaTop), c(
        "PCA for 10 samples using top 10 most variable windows.",
        "Using beta values and all 10 samples to calculate std dev."
    ))

    pcaSub <- getPCA(exampleTumourNormal, topVarNum = 10,
        topVarSamples = "_T", verbose = FALSE)
    expect_equal(plotLabels(pcaSub)[2],
        "Using beta values and 5 samples to calculate std dev.")

    grouped <- exampleTumourNormal %>%
        mutate(group = stringr::str_remove(sample_name, "_[NT]$"))
    pcaGroup <- getPCA(grouped, useGroupMeans = TRUE, topVarNum = 10,
        verbose = FALSE)
    expect_output(show(pcaGroup),
        "1 dimensionality reduction objects for 5 samples")
    expect_equal(plotLabels(pcaGroup), c(
        "PCA for 5 samples using top 10 most variable windows.",
        "Using beta values and all 5 samples to calculate std dev."
    ))

    set.seed(1)
    umap <- getUMAP(exampleTumourNormal, n_neighbors = 5, verbose = FALSE)
    expect_output(show(umap),
        "1 dimensionality reduction objects for 10 samples")
    expect_output(show(getResults(umap)$umap1),
        "UMAP result for 10 samples calculated over 599 windows")
    expect_equal(plotLabels(umap)[1],
        "UMAP for 10 samples using all 599 windows.")
})

test_that("getShapeScale() picks default shapes by number of levels", {

    shapes <- function(n, withNA = FALSE) {
        x <- tibble::tibble(s = letters[seq_len(n)])
        if (withNA) x$s[1] <- NA
        x
    }
    # the manual palette returns all its values: keep the first n
    palette <- function(scale, n) head(unname(scale$palette(n)), n)

    # NULLshape just wraps the palette and hides the legend
    nullScale <- getShapeScale(tibble::tibble(), "NULLshape", shapePalette = 21)
    expect_equal(palette(nullScale, 1), 21)
    expect_equal(nullScale$guide, "none")

    # few levels, no NA: filled shapes with a border
    scale <- getShapeScale(shapes(3), "s")
    expect_equal(palette(scale, 3), c(21, 24, 22))
    expect_true(is.na(scale$na.value))

    # NAs leave room for the NA shape (25)
    scale <- getShapeScale(shapes(3, withNA = TRUE), "s")
    expect_equal(palette(scale, 2), c(21, 24))
    expect_equal(scale$na.value, 25)

    # more levels fall back to unfilled shapes; the NA shape is then 7
    scale <- getShapeScale(shapes(7), "s")
    expect_equal(palette(scale, 4), c(16, 4, 0, 17))
    scale <- getShapeScale(shapes(6, withNA = TRUE), "s")
    expect_equal(palette(scale, 3), c(16, 4, 0))
    expect_equal(scale$na.value, 7)

    expect_error(
        getShapeScale(shapes(20), "s"),
        paste0("has 20 unique values; without NAs present, the maximum ",
            "allowed by default is 19")
    )
    expect_error(
        getShapeScale(shapes(20, withNA = TRUE), "s"),
        "with NAs present, the maximum allowed by default is 18"
    )
    # a diverging colour scale limits the default to filled shapes
    expect_error(
        getShapeScale(shapes(6), "s", colourScaleType = "diverging"),
        "maximum allowed by default when using a divergent colour scale is 5"
    )
})

test_that("getShapeScale() checks a custom shapePalette", {

    x <- tibble::tibble(s = c("a", "b", "c"))
    # the manual palette returns all its values: keep the first n
    palette <- function(scale, n) head(unname(scale$palette(n)), n)

    expect_equal(palette(getShapeScale(x, "s", c(1, 2, 3)), 3), c(1, 2, 3))
    expect_error(
        getShapeScale(x, "s", c(1, 2)),
        "`shapePalette` only has 2 values"
    )
    expect_error(
        getShapeScale(x, "s", c(1, 21, 22)),
        "either 0-20 \\(unfilled\\) OR 21-25 \\(filled\\), not both"
    )
    expect_error(
        getShapeScale(x, "s", list(1, 2, 3)),
        "must be 'line-first', 'filled-first', 'mixture', 'filled\\+border'"
    )

    expect_equal(
        palette(getShapeScale(x, "s", "filled+border", "qualitative"), 3),
        c(21, 24, 22)
    )
    expect_equal(
        palette(getShapeScale(x, "s", "line-first", "qualitative"), 3),
        c(1, 8, 2)
    )
    expect_equal(
        palette(getShapeScale(x, "s", "filled-first", "qualitative"), 3),
        c(16, 17, 15)
    )
    expect_equal(
        palette(getShapeScale(x, "s", "mixture", "qualitative"), 3),
        c(16, 4, 0)
    )
    expect_error(
        getShapeScale(x, "s", "bogus", "qualitative"),
        "must be 'line-first', 'filled-first', 'mixture', 'filled\\+border'"
    )
    expect_warning(
        getShapeScale(x, "s", "mixture", "diverging"),
        "Using 'mixture' palette with a divergent colour scale"
    )
    big <- tibble::tibble(s = letters[1:19])
    expect_error(
        getShapeScale(big, "s", "line-first", "qualitative"),
        "only 18 shapes available with `shapePalette` = 'line-first'"
    )
})

test_that("getShapeScale() resolves clashes with the NA shape", {

    x <- tibble::tibble(s = c("a", "b", "c", NA))
    # the manual palette returns all its values: keep the first n
    palette <- function(scale, n) head(unname(scale$palette(n)), n)

    # an NA shape from the wrong (filled/unfilled) set is reset
    expect_warning(
        scale <- getShapeScale(x, "s", NAshape = 7),
        "using filled shapes, but specified NAshape = 7 is from the unfilled"
    )
    expect_equal(scale$na.value, 25)

    expect_warning(
        scale <- getShapeScale(x, "s", c(1, 2, 3), NAshape = 21),
        "Resetting NAshape to default = 7 for unfilled shapes"
    )
    expect_equal(scale$na.value, 7)

    # an NA shape already in the palette is swapped out of it
    expect_warning(
        scale <- getShapeScale(x, "s", NAshape = 24),
        "Replacing that entry in the shapePalette"
    )
    expect_equal(scale$na.value, 24)
    expect_equal(palette(scale, 3), c(21, 22, 23))

    # fixed-order palettes reset the NA shape instead
    expect_warning(
        scale <- getShapeScale(
            x, "s", "line-first", "qualitative", NAshape = 1
        ),
        "has a fixed order, NAshape has been reset to default = 7"
    )
    expect_equal(scale$na.value, 7)

    # the default NA shape cannot also be a plotted shape
    expect_error(
        getShapeScale(x, "s", c(25, 21, 22)),
        "The default NAshape = 25 is already in use"
    )
    expect_error(
        suppressWarnings(getShapeScale(x, "s", c(7, 1, 2), NAshape = 21)),
        "The default NAshape = 7 is also already in use"
    )
})

test_that("getColourScale() builds the scale for each colour type", {

    filled <- ggplot2::scale_shape_manual(values = c(21, 24))
    unfilled <- ggplot2::scale_shape_manual(values = c(1, 2))
    x <- tibble::tibble(
        grp = c("a", "b", "c"), pos = c(1, 2, 3), neg = c(-3, -2, -1),
        div = c(-1, 0, 3)
    )
    scaleOf <- function(cV, type, shape = filled, cols = NULL, sym = FALSE) {
        getColourScale(x, cV, cols, type, shape, "grey50", sym)
    }

    expect_equal(scaleOf("NULLcol", "qualitative")$aesthetics, "fill")
    expect_equal(scaleOf("NULLcol", "qualitative", unfilled)$aesthetics,
        "colour")

    scale <- scaleOf("grp", "qualitative")
    expect_equal(scale$aesthetics, "fill")
    expect_equal(scale$na.value, "grey50")
    expect_equal(scaleOf("grp", "qualitative", unfilled)$aesthetics, "colour")

    for (type in c("sequential_non_neg", "sequential_non_pos")) {
        expect_s3_class(scaleOf("pos", type), "ScaleContinuous")
        expect_equal(scaleOf("pos", type, unfilled)$aesthetics, "colour")
    }
    # non-positive data runs the viridis scale in the other direction
    expect_false(identical(
        scaleOf("neg", "sequential_non_pos")$palette(c(0, 1)),
        scaleOf("pos", "sequential_non_neg")$palette(c(0, 1))
    ))

    # custom palettes
    scale <- scaleOf("grp", "qualitative", cols = c("red", "green", "blue"))
    expect_equal(unname(scale$palette(3)), c("red", "green", "blue"))
    scale <- scaleOf("grp", "qualitative", unfilled, c("red", "green", "blue"))
    expect_equal(scale$aesthetics, "colour")
    expect_error(
        scaleOf("grp", "qualitative", cols = c("red", "green")),
        "`colour` variable 'grp' has 3 unique values; the `colourPalette`"
    )
    expect_s3_class(
        scaleOf("pos", "sequential_non_neg", cols = c("white", "black")),
        "ScaleContinuous"
    )
    expect_equal(
        scaleOf("pos", "sequential_non_pos", unfilled,
            c("white", "black"))$aesthetics,
        "colour"
    )
})

test_that("getColourScale() centres diverging scales on zero", {

    filled <- ggplot2::scale_shape_manual(values = c(21, 24))
    unfilled <- ggplot2::scale_shape_manual(values = c(1, 2))
    scaleOf <- function(values, shape = filled, sym = FALSE) {
        getColourScale(
            tibble::tibble(div = values), "div", NULL, "diverging", shape,
            "grey50", sym
        )
    }

    for (sym in c(FALSE, TRUE)) {
        # larger positive and larger negative extremes take different branches
        for (values in list(c(-1, 0, 3), c(-3, 0, 1))) {
            scale <- scaleOf(values, sym = sym)
            expect_s3_class(scale, "ScaleContinuous")
            expect_equal(scale$aesthetics, "fill")
            expect_equal(scale$na.value, "grey50")
        }
        expect_equal(scaleOf(c(-1, 0, 3), unfilled, sym)$aesthetics, "colour")
    }

    # grey90 sits at the midpoint of the RdBu palette: zero maps to it
    scale <- scaleOf(c(-2, 0, 2), sym = TRUE)
    expect_match(scale$palette(0.5), "^#E5E5E5", ignore.case = TRUE)
})

test_that("getPCA() rejects bad topVarNum, topVarSamples and methods", {

    pca <- function(topVarNum = 50, ...) {
        getPCA(exampleTumourNormal, topVarNum = topVarNum, verbose = FALSE, ...)
    }

    expect_error(
        pca(topVarNum = "a"), "topVarNum should be a numeric vector"
    )
    expect_error(
        pca(topVarSamples = 1:3),
        "topVarSamples should be NULL, a character vector"
    )
    expect_error(
        pca(topVarSamples = c("Colon1_N", "bogus")),
        "topVarSamples contains sample names that are not in the qseaSet"
    )
    expect_error(
        pca(topVarNum = c(50, 100), topVarSamples = list("_T$", "_N$", "x")),
        "If topVarSamples is a list and length\\(topVarNum\\) > 1"
    )
    expect_error(
        getDimRed(exampleTumourNormal, method = "foo", verbose = FALSE),
        "Method foo not known! Options are PCA or UMAP."
    )
})

test_that("getPCA() takes regions, nrpm values and per-set topVarSamples", {

    pca <- function(...) {
        getPCA(exampleTumourNormal, verbose = FALSE, ...)
    }
    regions <- qsea::getRegions(exampleTumourNormal)[1:300]

    within <- pca(regionsToOverlap = regions, topVarNum = 50)
    expect_equal(dim(getCoordinates(within)$pca1), c(10, 5))
    expect_equal(getParameters(within)$topVar$topVarNum, 50)

    nrpm <- pca(normMethod = "nrpm", topVarNum = 50)
    expect_equal(getParameters(nrpm)$normMethod, "nrpm")
    expect_false(isTRUE(all.equal(
        getCoordinates(nrpm)$pca1,
        getCoordinates(pca(topVarNum = 50))$pca1
    )))

    # a list gives one PCA per topVarNum, each on its own samples
    perSet <- pca(topVarNum = c(50, 100), topVarSamples = list("_T$", "_N$"))
    expect_length(getResults(perSet), 2)

    scaled <- pca(topVarNum = 50, center = FALSE, scale = TRUE, nPC = 3)
    expect_equal(dim(getCoordinates(scaled)$pca1), c(10, 3))
})

test_that("getPCA() checks a supplied dataTable against the qseaSet", {

    betas <- getBetaTable(exampleTumourNormal)

    renamed <- betas
    colnames(renamed)[colnames(renamed) == "Colon1_T"] <- "nope"
    expect_error(
        getPCA(exampleTumourNormal, dataTable = renamed, topVarNum = 50,
            verbose = FALSE),
        "does not have a matching"
    )

    # values that differ from what the qseaSet gives for the same windows
    shuffled <- betas
    shuffled$Colon1_T <- rev(shuffled$Colon1_T)
    expect_warning(
        getPCA(exampleTumourNormal, dataTable = shuffled, topVarNum = 50,
            verbose = FALSE),
        "does not match input dataTable"
    )
})

test_that("plotDimRed() validates colour, shape and palettes", {

    obj <- getPCA(exampleTumourNormal, topVarNum = 100, verbose = FALSE)

    expect_error(plotDimRed("a"), "First argument should be the output from")
    expect_error(
        plotDimRed(obj, colour = "nope"),
        "Can't colour by nope, as it is not present in the sampleTable!"
    )
    expect_error(
        plotDimRed(obj, shape = "nope"),
        "Can't set shapes by nope, as it is not present in the sampleTable!"
    )
    expect_error(
        plotDimRed(obj, colourPalette = c("red", "blue")),
        "`colourPalette` argument is non-NULL, but `colour` argument is NULL."
    )
    expect_error(
        plotDimRed(obj, shape = c("type", "tissue")),
        "Argument `shape` can only be of length one."
    )

    logicalObj <- obj
    logicalObj@sampleTable$isTumour <- logicalObj@sampleTable$tumour == "Tumour"
    expect_error(
        plotDimRed(logicalObj, colour = "isTumour"),
        "The variable `isTumour` can not be mapped to a colour scale."
    )
})

test_that("plotDimRed() builds one plot per colour variable and component", {

    obj <- getPCA(exampleTumourNormal, topVarNum = 100, verbose = FALSE)

    plots <- plotDimRed(obj, colour = c("tissue", "age"), components = c(1, 3))
    expect_named(plots, "pca1")
    expect_named(plots$pca1, c("tissue", "age"))
    expect_named(plots$pca1$tissue, "PC1vsPC3")
    expect_s3_class(plots$pca1$age$PC1vsPC3, "ggplot")

    # no colour: the plots sit directly under the result
    plain <- plotDimRed(obj)
    expect_named(plain$pca1, c("PC1vsPC2", "PC2vsPC3"))

    labelled <- plotDimRed(
        obj, colour = "type", shape = "tissue", showSampleNames = TRUE,
        colourPalette = c("red", "blue", "green", "orange", "purple"),
        shapePalette = "filled+border", symDivColourScale = TRUE
    )
    expect_s3_class(labelled$pca1$type$PC1vsPC2, "ggplot")

    # without ggrepel, the labels fall back to geom_text()
    local_mocked_bindings(requireNamespace = function(...) FALSE,
        .package = "base")
    expect_message(
        plotDimRed(obj, showSampleNames = TRUE),
        "Package \"ggrepel\" is recommended"
    )
})
