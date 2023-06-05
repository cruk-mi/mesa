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
    
    # Test that the x matrices in the res components are not empty
    expect_false(is_empty(obj1$res$pca1$x))
    expect_false(is_empty(obj2$res$pca1$x))
    expect_false(is_empty(obj3$res$pca1$x))
    expect_false(is_empty(obj4$res$pca1$x))
    expect_false(is_empty(obj5$res$pca1$x))
    expect_false(is_empty(obj6$res$pca1$x))
    expect_false(is_empty(obj7$res$pca1$x))
    
    # Test that the dimensions of the x matrices are consistent
    expect_equal(dim(obj1$res$pca1$x), dim(obj2$res$pca1$x))
    expect_equal(dim(obj1$res$pca1$x), dim(obj3$res$pca1$x))
    expect_equal(dim(obj1$res$pca1$x), dim(obj4$res$pca1$x))
    expect_equal(dim(obj1$res$pca1$x), dim(obj5$res$pca1$x))
    expect_equal(dim(obj1$res$pca1$x), dim(obj6$res$pca1$x))
    expect_equal(dim(obj1$res$pca1$x), dim(obj7$res$pca1$x))
    
    # Test that the variance_explained component is a numeric vector
    expect_true(is.numeric(obj1$res$pca1$variance_explained))
    
    # Test that getPCA() returns an error when called with invalid arguments
    expect_error(exampleTumourNormal %>% getPCA("invalid_argument"))
    
    # Check res components are of correct length
    expect_equal(length(obj1$res), 1)
    expect_equal(length(obj6$res), 3)
    
    # Check PCA results are different for different numbers of windows
    expect_false(isTRUE(all.equal(obj6$res$pca1$x,obj6$res$pca2$x)))
    expect_false(isTRUE(all.equal(obj6$res$pca1$x,obj6$res$pca3$x)))
    
    # Check dimensions of x are correct
    expect_equal(dim(obj6$res$pca1$x), c(10, 5))
    expect_equal(dim(obj6$res$pca2$x), c(10, 5))
    expect_equal(dim(obj6$res$pca3$x), c(10, 5))
    
    # Check dimensions of rotation are correct
    expect_equal(dim(obj6$res$pca1$rotation), c(10, 5))
    expect_equal(dim(obj6$res$pca2$rotation), c(100, 5))
    expect_equal(dim(obj6$res$pca3$rotation), c(200, 5))
    
    ##### plotting tests
    
    # Check no error is returned
    expect_no_error(plots1 <- plotPCA(obj1, exampleTumourNormal))
    expect_no_error(plots6 <- plotPCA(obj6, exampleTumourNormal))
    
    # Check output is correct length
    expect_equal(length(plots1), 1)
    expect_equal(length(plots6), 3)
    
    # Test that plotPCA() returns a ggplot object
    expect_equal(class(plots1 <- plotPCA(obj1, exampleTumourNormal)), "ggplot")
    expect_equal(class(plots6 <- plotPCA(obj6, exampleTumourNormal)), "ggplot")
    
    # Test other arguments of plotPCA()
    expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "type"))
    expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "age"))
    expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "gender"))
    expect_no_error(plotPCA(obj1, exampleTumourNormal %>% mutate(diver = seq(-4,5)), colour = "diver"))
    expect_no_error(plotPCA(obj1, exampleTumourNormal, shape = "type"))
    
    expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "type", colourPalette = RColorBrewer::brewer.pal(5,"Oranges")))
    expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "gender", shapePalette = c(2,5), shape = "gender"))


  # randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000000) %>%
  #   mutate(patient = stringr::str_remove(sample_name,"[TN]$"),
  #          variableWithOneLevel = "Test",
  #          experiment = ifelse( stringr::str_detect(sample_name,"[1234]"),"A","B"),
  #          experimentConfounded = ifelse( stringr::str_detect(sample_name,"[1234]N"),"A","B"),
  #          fish = c(9,1,2,5,3,4,1,5,-1,2,-1,-1,-2,0,4,3)
  #   )

}

)

test_that("UMAPs", {

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

  expect_true(isTRUE(all.equal(obj1$res$umap1$x,obj2$res$umap1$x)))
  expect_false(isTRUE(all.equal(obj3$res$umap1$x,obj4$res$umap1$x)))
  expect_false(isTRUE(all.equal(obj1$res$umap1$x,obj5$res$umap1$x)))
  expect_false(isTRUE(all.equal(obj1$res$umap1$x,obj6$res$umap1$x)))
  expect_true(isTRUE(all.equal(obj1$res$umap1$x,obj7$res$umap1$x)))

  expect_equal(length(obj1$res), 1)

  expect_equal(length(obj6$res), 3)
  expect_false(isTRUE(all.equal(obj6$res$umap1$x,obj6$res$umap2$x)))
  expect_false(isTRUE(all.equal(obj6$res$umap1$x,obj6$res$umap3$x)))

  expect_no_error(plots1 <- plotUMAP(obj1, exampleTumourNormal))
  expect_no_error(plots6 <- plotUMAP(obj6, exampleTumourNormal))

  expect_no_error(plotUMAP(obj1, exampleTumourNormal, colour = "type"))
  expect_no_error(plotUMAP(obj1, exampleTumourNormal, colour = "age"))
  expect_no_error(plotUMAP(obj1, exampleTumourNormal, colour = "gender"))
  expect_no_error(plotUMAP(obj1, exampleTumourNormal %>% mutate(diver = seq(-4,5)), colour = "diver"))
  expect_no_error(plotUMAP(obj1, exampleTumourNormal, shape = "type"))

  expect_no_error(plotUMAP(obj1, exampleTumourNormal, colour = "type", colourPalette = RColorBrewer::brewer.pal(5,"Oranges")))
  expect_no_error(plotUMAP(obj1, exampleTumourNormal, colour = "gender", shapePalette = c(2,5), shape = "gender"))

  expect_equal(length(plots1), 1)
  expect_equal(length(plots6), 3)

  # randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000000) %>%
  #   mutate(patient = stringr::str_remove(sample_name,"[TN]$"),
  #          variableWithOneLevel = "Test",
  #          experiment = ifelse( stringr::str_detect(sample_name,"[1234]"),"A","B"),
  #          experimentConfounded = ifelse( stringr::str_detect(sample_name,"[1234]N"),"A","B"),
  #          fish = c(9,1,2,5,3,4,1,5,-1,2,-1,-1,-2,0,4,3)
  #   )

}

)
