test_that("Parallelisation works with useParallel", {
  
  expect_false(getMesaParallel())
  
  expect_true(setMesaParallel(useParallel = TRUE))
  
  expect_true(getMesaParallel())
  
  expect_false(setMesaParallel(useParallel = FALSE))

  expect_false(getMesaParallel())  
}
)

test_that("Parallelisation works with nCores", {
  
  expect_false(getMesaParallel())
  
  expect_true(setMesaParallel(nCores = 2))
  
  expect_true(getMesaParallel())
  
}
)
