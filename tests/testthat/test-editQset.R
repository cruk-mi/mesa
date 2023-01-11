test_that("Filtering Qset works", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 1000)

  expect_equal(randomSet %>%
                 filter(stringr::str_detect(sample_name, "Sim7")) %>%
                 qsea::getSampleNames(), c("Sim7T","Sim7N"))

  expect_equal(randomSet %>%
                 filter(group == "Tumor") %>%
                 qsea::getSampleNames() %>%
                 length(), 8)
  })

test_that("Reducing Qset Windows", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000)

  regions <- qsea::getRegions(randomSet)[1:1000]

  expect_equal(randomSet %>%
                 filterByOverlaps(regions) %>%
                 qsea::getRegions() %>% length(), 1000)

})


test_that("Mutating Qset works", {

  sampTab <- qsea::getExampleQseaSet(repl = 3, expSamplingDepth = 1000) %>%
    mutate(newCol = ifelse(stringr::str_detect(sample_name,"Sim1"),"fish","duck")) %>%
    mutate(repNum = stringr::str_remove(sample_name,"Sim")) %>%
    qsea::getSampleTable()

  expect_equal(
    sampTab %>% pull(newCol) %>% unique(), c("fish","duck")
  )

  expect_equal(
    sampTab %>% pull(repNum), c("1T","2T","3T","1N","2N","3N")
  )

})

test_that("Sorting Qsets", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000)

  sortedSet <- randomSet %>%
    sort()

  expect_equal(randomSet %>% getSampleNames() %>% sort() , sortedSet %>% getSampleNames())

})

test_that("Combining Qsets", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000)

  splitSet <- filter(randomSet, group == "Tumor") %>%
    combineQsets(filter(randomSet, group != "Tumor"))

  expect_equal(splitSet %>% sort(), randomSet %>% sort() )

})

test_that("filterByOverlaps", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000)

  reducedGRanges <- randomSet %>% qsea::getRegions() %>% filter(start <= 1000000)

  expect_equal(randomSet %>%
    filterByOverlaps(reducedGRanges) %>%
    qsea::getRegions() %>%
    length(), 2000)

  expect_equal(randomSet %>%
                 filterByNonOverlaps(reducedGRanges) %>%
                 qsea::getRegions() %>%
                 length(), 8000)

})

test_that("addLibraryInformation", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000)

  colNames <- randomSet %>%
    addLibraryInformation() %>%
    qsea::getSampleTable() %>%
    colnames()

  expect_true("fragment_length" %in% colNames)
  expect_true("valid_fragments" %in% colNames)

})


test_that("renameSamples", {

  randomSet <- qsea::getExampleQseaSet(repl = 8, expSamplingDepth = 100000) %>%
    mutate(new_column = LETTERS[1:16])

  renamedSet <- randomSet %>%
    renameSamples("new_column")

  expect_true(all(getSampleNames(renamedSet)  ==  LETTERS[1:16]))

})

test_that("is.qseaSet", {

  expect_true(qsea::getExampleQseaSet() %>%
                is.qseaSet())

  expect_false(qsea::getExampleQseaSet() %>%
                 qsea::getSampleTable() %>%
                 is.qseaSet())

  expect_false(is.qseaSet("string"))
  expect_false(is.qseaSet(2))
  expect_false(is.qseaSet(1:10))

})
