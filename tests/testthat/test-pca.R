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
    expect_no_error(print(plotPCA(obj1 %>% mutate(PC1 = 0))))
    expect_no_error(plotPCA(obj1, shape = "type"))

    expect_no_error(plotPCA(obj1, colour = "type", colourPalette = RColorBrewer::brewer.pal(5,"Oranges")))
    expect_no_error(plotPCA(obj1, colour = "gender", shapePalette = c(2,5), shape = "gender"))

    expect_no_error(plotPCA(obj1 %>% mutate(gender = as.factor(gender)),
                            colour = "gender", shapePalette = c(2,5), shape = "gender"))

    expect_no_error(obj1 %>% mutate(newCol = rnorm(10)) %>% plotPCA(colour = "newCol"))

}

)

test_that("UMAPs", {

  skip_long_checks()

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
