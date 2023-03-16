test_that("PCAs", {

  expect_no_error(obj1 <- exampleTumourNormal %>% getPCA() )
  expect_no_error(obj2 <- exampleTumourNormal %>% getPCA(returnDataTable = TRUE) )
  expect_no_error(obj3 <- exampleTumourNormal %>% getPCA(topVarNum = 10) )
  expect_no_error(obj4 <- exampleTumourNormal %>% getPCA(topVarSamples = "_T", topVarNum = 10) )
  expect_no_error(obj5 <- exampleTumourNormal %>% getPCA(minDensity = 10) )
  expect_no_error(obj6 <- exampleTumourNormal %>% getPCA(topVarNum = c(10,100,200)) )


  expect_no_error(obj7 <- exampleTumourNormal %>% getPCA(dataTable = getBetaTable(exampleTumourNormal)) )

  expect_true(isTRUE(all.equal(obj1$pca$pca1$x,obj2$pca$pca1$x)))
  expect_false(isTRUE(all.equal(obj1$pca$pca1$x,obj3$pca$pca1$x)))
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

  expect_no_error(plotPCA(obj1, exampleTumourNormal, colourVar = "type"))
  expect_no_error(plotPCA(obj1, exampleTumourNormal, colourVar = "age"))
  expect_no_error(plotPCA(obj1, exampleTumourNormal, colourVar = "gender"))
  expect_no_error(plotPCA(obj1, exampleTumourNormal %>% mutate(diver = seq(-4,5)), colourVar = "diver"))
  expect_no_error(plotPCA(obj1, exampleTumourNormal, shapeVar = "type"))

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

exampleTumourNormal