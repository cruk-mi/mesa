test_that("PCAs", {

  expect_no_error(obj1 <- exampleTumourNormal %>% getPCA() )
  expect_no_error(obj2 <- exampleTumourNormal %>% getPCA(returnDataTable = TRUE) )
  expect_no_error(obj3 <- exampleTumourNormal %>% getPCA(topVarNum = 10) )
  expect_no_error(obj4 <- exampleTumourNormal %>% getPCA(topVarSamples = "_T", topVarNum = 10) )
  expect_no_error(obj5 <- exampleTumourNormal %>% getPCA(minDensity = 10) )
  expect_no_error(obj6 <- exampleTumourNormal %>% getPCA(topVarNum = c(10,100,200)) )


  expect_no_error(obj7 <- exampleTumourNormal %>% getPCA(dataTable = getBetaTable(exampleTumourNormal)) )

  expect_true(isTRUE(all.equal(obj1$pca$pca1$x,obj2$pca$pca1$x)))
  expect_false(isTRUE(all.equal(obj3$pca$pca1$x,obj4$pca$pca1$x)))
  expect_false(isTRUE(all.equal(obj1$pca$pca1$x,obj5$pca$pca1$x)))
  expect_false(isTRUE(all.equal(obj1$pca$pca1$x,obj6$pca$pca1$x)))
  expect_true(isTRUE(all.equal(obj1$pca$pca1$x,obj7$pca$pca1$x)))

  expect_equal(length(obj1$pca), 1)

  expect_equal(length(obj6$pca), 3)
  expect_false(isTRUE(all.equal(obj6$pca$pca1$x,obj6$pca$pca2$x)))
  expect_false(isTRUE(all.equal(obj6$pca$pca1$x,obj6$pca$pca3$x)))

  expect_no_error(plots1 <- plotPCA(obj1, exampleTumourNormal))
  expect_no_error(plots6 <- plotPCA(obj6, exampleTumourNormal))

  expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "type"))
  expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "age"))
  expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "gender"))
  expect_no_error(plotPCA(obj1, exampleTumourNormal %>% mutate(diver = seq(-4,5)), colour = "diver"))
  expect_no_error(plotPCA(obj1, exampleTumourNormal, shape = "type"))

  expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "type", colourPalette = RColorBrewer::brewer.pal(5,"Oranges")))
  expect_no_error(plotPCA(obj1, exampleTumourNormal, colour = "gender", shapePalette = c(2,5), shape = "gender"))

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

test_that("UMAPs", {

  expect_no_error(obj1 <- exampleTumourNormal %>% getUMAP() )
  expect_no_error(obj2 <- exampleTumourNormal %>% getUMAP(returnDataTable = TRUE) )
  expect_no_error(obj3 <- exampleTumourNormal %>% getUMAP(topVarNum = 10) )
  expect_no_error(obj4 <- exampleTumourNormal %>% getUMAP(topVarSamples = "_T", topVarNum = 10) )
  expect_no_error(obj5 <- exampleTumourNormal %>% getUMAP(minDensity = 10) )
  expect_no_error(obj6 <- exampleTumourNormal %>% getUMAP(topVarNum = c(10,100,200)) )


  expect_no_error(obj7 <- exampleTumourNormal %>% getUMAP(dataTable = getBetaTable(exampleTumourNormal)) )

  expect_true(isTRUE(all.equal(obj1$umap$umap1$x,obj2$umap$umap1$x)))
  expect_false(isTRUE(all.equal(obj3$umap$umap1$x,obj4$umap$umap1$x)))
  expect_false(isTRUE(all.equal(obj1$umap$umap1$x,obj5$umap$umap1$x)))
  expect_false(isTRUE(all.equal(obj1$umap$umap1$x,obj6$umap$umap1$x)))
  expect_true(isTRUE(all.equal(obj1$umap$umap1$x,obj7$umap$umap1$x)))

  expect_equal(length(obj1$umap), 1)

  expect_equal(length(obj6$umap), 3)
  expect_false(isTRUE(all.equal(obj6$umap$umap1$x,obj6$umap$umap2$x)))
  expect_false(isTRUE(all.equal(obj6$umap$umap1$x,obj6$umap$umap3$x)))

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