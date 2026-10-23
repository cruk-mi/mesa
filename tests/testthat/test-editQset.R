test_that("Filtering Qset works", {

    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 1000)

    expect_equal(randomSet %>%
        filter(stringr::str_detect(sample_name, "Sim7")) %>%
        qsea::getSampleNames(), c("Sim7T", "Sim7N"))

    expect_equal(randomSet %>%
        filter(group == "Tumor") %>%
        qsea::getSampleNames() %>%
        length(), 8)
})

test_that("Changing Qset Windows", {

    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 100000)

    regions <- qsea::getRegions(randomSet)[1:1000]

    expect_equal(randomSet %>%
        filterByOverlaps(regions) %>%
        qsea::getRegions() %>% length(), 1000)

    expect_equal(randomSet %>%
        filterByNonOverlaps(regions) %>%
        qsea::getRegions() %>% length(), 9000)
})


test_that("Mutating Qset works", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    sampTab <- qseaSet %>%
        mutate(newCol = ifelse(stringr::str_detect(sample_name, "Sim1"), "fish", "duck")) %>%
        mutate(repNum = stringr::str_remove(sample_name, "Sim")) %>%
        qsea::getSampleTable()

    expect_equal(
        sampTab %>% pull(newCol) %>% unique(), c("fish", "duck")
    )

    expect_equal(
        sampTab %>% pull(repNum), c("1T", "2T", "3T", "1N", "2N", "3N")
    )

    expect_error(
        qseaSet %>% mutate(sample_name = "new_name")
    )

    expect_error(
        qseaSet %>% mutate(sample_name = "new_name")
    )

})


test_that("Join works", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    expect_true(
        "new" %in% (colnames(qseaSet %>% left_join(tibble(group = c("Tumor", "Normal"), new = 1:2)) %>% getSampleTable()))
    )

    expect_true(
        "new" %in% (colnames(qseaSet %>% left_join(tibble(not_group = c("Tumor", "Normal"),
            new = 1:2),
        by = c("group" = "not_group")) %>% getSampleTable()))
    )

})


test_that("Sorting Qsets", {

    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 100000)

    sortedSet <- randomSet %>%
        sort()

    expect_equal(randomSet %>% getSampleNames() %>% sort(), sortedSet %>% getSampleNames())

})

test_that("Combining Qsets", {

    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 1000)

    splitSet <- filter(randomSet, group == "Tumor") %>%
        combineQsets(filter(randomSet, group != "Tumor"))

    splitSetList <- list(a = filter(randomSet, group == "Tumor"),
        b = filter(randomSet, group != "Tumor")) %>%
        combineQsetsList()

    expect_equal(splitSet %>% sort(), randomSet %>% sort())
    expect_equal(splitSetList %>% sort(), randomSet %>% sort())
    expect_equal(splitSet %>% sort(), splitSetList %>% sort())

})

test_that("filterByOverlaps", {

    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 100000)

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

    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 100000)

    colNames <- randomSet %>%
        addLibraryInformation() %>%
        qsea::getSampleTable() %>%
        colnames()

    expect_true("fragment_length" %in% colNames)
    expect_true("valid_fragments" %in% colNames)

})


test_that("renameSamples", {

    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 1000) %>%
        mutate(new_column = LETTERS[1:16])

    renamedSet <- randomSet %>%
        renameSamples("new_column")

    expect_true(all(getSampleNames(renamedSet) == LETTERS[1:16]))

    expect_true("Example1T" %in% (randomSet %>% renameQsetNames("Sim", "Example") %>% getSampleNames()))

})

test_that("poolSamples works on the example data (#152)", {

    # The pre-filtered example windows are all above CpG density 5, so
    # addNormalisation() has no background windows to fit.
    pooled <- poolSamples(
        exampleTumourNormal, "_[TN]$", renormalise = FALSE
    )
    patients <- c("Colon1", "Colon2", "Lung1", "Lung2", "Lung3")

    expect_equal(getSampleNames(pooled), patients)
    expect_equal(
        unname(qsea::getCounts(pooled)[, "Colon1"]),
        unname(rowSums(
            qsea::getCounts(exampleTumourNormal)[, c("Colon1_N", "Colon1_T")]
        ))
    )
    expect_equal(rownames(pooled@enrichment$parameters), patients)
    expect_equal(colnames(pooled@enrichment$factors), patients)
    expect_equal(rownames(pooled@zygosity), patients)

    # every library column survives pooling, not just the matched ones
    expect_setequal(
        colnames(pooled@libraries$file_name),
        colnames(exampleTumourNormal@libraries$file_name)
    )
    expect_true(all(is.na(pooled@libraries$file_name$offset)))
    expect_equal(nrow(getSampleQCSummary(pooled)), 5L)
})

test_that("poolSamples works on a set without input libraries (#152)", {

    # qsea's example set has no input_file library table
    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 100000)
    expect_null(randomSet@libraries$input_file)

    pooled <- poolSamples(randomSet, "[TN]$", maxPatternDensity = 0.5)

    expect_equal(getSampleNames(pooled), paste0("Sim", 1:8))
    expect_equal(dim(pooled@zygosity), c(8L, 1L))
    expect_false(anyNA(pooled@libraries$file_name$offset))
    expect_true("fragment_length" %in% colnames(pooled@libraries$file_name))
    expect_equal(ncol(getBetaTable(pooled)), 4L + 8L)
})

test_that("poolSamples pools only the samples that share a name (#152)", {

    # "Sim1" is a prefix of "Sim12", but the two must stay separate pools
    randomSet <- cachedExampleQset(repl = 8, expSamplingDepth = 100000) %>%
        renameQsetNames("^Sim2", "Sim12")
    pooled <- poolSamples(randomSet, "[TN]$", renormalise = FALSE)
    counts <- qsea::getCounts(randomSet)

    expect_true(all(c("Sim1", "Sim12") %in% getSampleNames(pooled)))
    expect_equal(
        unname(qsea::getCounts(pooled)[, "Sim1"]),
        unname(rowSums(counts[, c("Sim1T", "Sim1N")]))
    )
    expect_equal(
        unname(qsea::getCounts(pooled)[, "Sim12"]),
        unname(rowSums(counts[, c("Sim12T", "Sim12N")]))
    )
})

test_that("poolSamples keeps the enrichment method of the set (#152)", {

    noEnrichment <- cachedExampleQset(repl = 8, expSamplingDepth = 100000) %>%
        addNormalisation(enrichmentMethod = "none", maxPatternDensity = 0.5)
    pooled <- poolSamples(noEnrichment, "[TN]$", maxPatternDensity = 0.5)

    expect_equal(qsea::getParameters(pooled)$enrichmentMethod, "none")
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

test_that("select.qseaSet", {
    expect_equal(exampleTumourNormal %>% select(tissue) %>% getSampleTable() %>% colnames(), c("sample_name", "group", "tissue"))
    expect_equal(exampleTumourNormal %>% select(type) %>% getSampleTable() %>% colnames(), c("sample_name", "group", "type"))
    expect_contains(exampleTumourNormal %>% select(-starts_with("t")) %>% getSampleTable() %>% colnames(), c("sample_name", "group"))

    expect_equal(exampleTumourNormal %>% select(-starts_with("t")) %>% getSampleTable() %>% ncol(), 6)
    expect_equal(exampleTumourNormal %>% select(-age) %>% getSampleTable() %>% ncol(), 8)
    expect_contains(exampleTumourNormal %>% select(-age) %>% getSampleTable() %>% colnames(), c("sample_name", "group"))

})
