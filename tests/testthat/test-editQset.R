context("Merging and filtering qsets")

test_that("Filtering Qset works", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000)

  expect_equal(randomSet %>%
                 filterQset(stringr::str_detect(sample_name, "Sim7")) %>%
                 qsea::getSampleNames(), c("Sim7T","Sim7N"))

  expect_equal(randomSet %>%
                 filterQset(group == "Tumor") %>%
                 qsea::getSampleNames() %>%
                 length(), 8)
  })

test_that("Dropping Pooled Control works", {

  randomSet <- qsea::getExampleQseaSet(repl = 3, expSamplingDepth = 1000) %>%
    renameQsetNames("Sim3N", "PooledControl")

  expect_true("PooledControl" %in% qsea::getSampleNames(randomSet))

  expect_false("PooledControl" %in% qsea::getSampleNames(dropPooledControl(randomSet)))

})

test_that("Reducing Qset Windows", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000) %>%
    renameQsetNames("Sim8N", "PooledControl")

  regions <- qsea::getRegions(randomSet)[1:1000]


  expect_equal(randomSet %>%
                 filterByOverlaps(regions) %>%
                 qsea::getRegions() %>% length(), 1000)

})


test_that("Mutating Qset works", {

  sampTab <- qsea::getExampleQseaSet(repl = 3, expSamplingDepth = 1000) %>%
    mutateQset(newCol = ifelse(stringr::str_detect(sample_name,"Sim1"),"fish","duck")) %>%
    mutateQset(repNum = stringr::str_remove(sample_name,"Sim")) %>% getSampleTable()

  expect_equal(
    sampTab %>% pull(newCol) %>% unique(), c("fish","duck")
  )

  expect_equal(
    sampTab %>% pull(repNum), c("1T","2T","3T","1N","2N","3N")
  )

})

test_that("Sorting Qsets", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000) %>%
    renameQsetNames("Sim8N", "PooledControl")

  sortedSet <- randomSet %>%
    sortQset()

  expect_equal(randomSet %>% getSampleNames() %>% sort() , sortedSet %>% getSampleNames())

})

test_that("Combining Qsets", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000) %>%
    renameQsetNames("Sim8N", "PooledControl")

  splitSet <- filterQset(randomSet, group == "Tumor") %>%
    combineQsets(filterQset(randomSet, group != "Tumor"))

  expect_equal(splitSet %>% sortQset(), randomSet %>% sortQset() )

})

test_that("filterByOverlaps", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000) %>%
    renameQsetNames("Sim8N", "PooledControl")

  reducedGRanges <- randomSet %>% qsea::getRegions() %>% filter(start <= 1000000)

  expect_equal(randomSet %>%
    filterByOverlaps(reducedGRanges) %>%
    getRegions() %>%
    length(), 2000)

  expect_equal(randomSet %>%
                 filterByNonOverlaps(reducedGRanges) %>%
                 getRegions() %>%
                 length(), 8000)

})

test_that("addLibraryInformation", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000) %>%
    renameQsetNames("Sim8N", "PooledControl")

  colNames <- randomSet %>%
    addLibraryInformation() %>%
    getSampleTable() %>%
    colnames()

  expect_true("fragment_length" %in% colNames)
  expect_true("valid_fragments" %in% colNames)

})
