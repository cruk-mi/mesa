test_that("getLibraries() returns the whole list or one table", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    expect_type(getLibraries(qseaSet), "list")
    expect_true("file_name" %in% names(getLibraries(qseaSet)))

    libraryTable <- getLibraries(qseaSet, "file_name")
    expect_equal(getLibraries(qseaSet)$file_name, libraryTable)
    expect_equal(rownames(libraryTable), qsea::getSampleNames(qseaSet))
    expect_true("valid_fragments" %in% colnames(libraryTable))

    expect_null(getLibraries(qseaSet, "notATable"))
})


test_that("setLibraries() replaces one table or the whole list", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    libraryTable <- getLibraries(qseaSet, "file_name")
    libraryTable[, "valid_fragments"] <- 1

    newSet <- setLibraries(qseaSet, libraryTable, "file_name")
    expect_equal(
        unname(getLibraries(newSet, "file_name")[, "valid_fragments"]),
        rep(1, length(qsea::getSampleNames(qseaSet)))
    )

    expectedLibraries <- getLibraries(qseaSet)
    expectedLibraries$file_name <- libraryTable
    expect_equal(getLibraries(newSet), expectedLibraries)

    newLibraries <- list(file_name = libraryTable)
    newSet <- setLibraries(qseaSet, newLibraries)
    expect_equal(getLibraries(newSet), newLibraries)

    inputTable <- getLibraries(qseaSet, "file_name")
    libraryTable[, "valid_fragments"] <- 2
    libraryList <- list("file_name" = libraryTable, "input_file" = inputTable)
    newSet <- setLibraries(qseaSet, libraryList)

    expect_equal(getLibraries(newSet), libraryList)

})


test_that("getEnrichment() returns the whole list or one element", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    expect_type(getEnrichment(qseaSet), "list")
    expect_equal(
        getEnrichment(qseaSet, "pattern_name"),
        getEnrichment(qseaSet)$pattern_name
    )
    expect_null(getEnrichment(qseaSet, "notAnElement"))
})


test_that("setEnrichment() replaces one element", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    newSet <- setEnrichment(qseaSet, "AT", "pattern_name")

    expect_equal(getEnrichment(newSet, "pattern_name"), "AT")
    expect_equal(
        getEnrichment(newSet, "factors"), getEnrichment(qseaSet, "factors")
    )
})


test_that("setCounts() and setRegions() replace the window data", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    counts <- qsea::getCounts(qseaSet)
    newSet <- setCounts(qseaSet, counts * 0)

    expect_equal(sum(qsea::getCounts(newSet)), 0)
    expect_equal(dim(qsea::getCounts(newSet)), dim(counts))

    regions <- qsea::getRegions(qseaSet)
    regions$CpG_density <- 1
    expect_equal(
        unique(qsea::getRegions(setRegions(qseaSet, regions))$CpG_density), 1
    )
})


test_that("setSampleTable() replaces the sample annotations", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    sampleTable <- qsea::getSampleTable(qseaSet)
    sampleTable$batch <- "A"

    newSet <- setSampleTable(qseaSet, sampleTable)

    expect_equal(
        qsea::getSampleTable(newSet)$batch,
        rep("A", length(qsea::getSampleNames(qseaSet)))
    )
    expect_equal(qsea::getSampleNames(newSet), qsea::getSampleNames(qseaSet))
})


test_that("setParameters() merges and overwrites by name", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)

    newSet <- setParameters(qseaSet, list(minMapQual = 10, newParam = "x"))

    expect_equal(qsea::getParameters(newSet)$minMapQual, 10)
    expect_equal(qsea::getParameters(newSet)$newParam, "x")

    overwrittenSet <- setParameters(newSet, list(newParam = "updated"))
    expect_equal(qsea::getParameters(overwrittenSet)$newParam, "updated")

    # pre-existing parameters survive
    expect_equal(
        qsea::getParameters(overwrittenSet)$window_size,
        qsea::getParameters(qseaSet)$window_size
    )

    expect_error(setParameters(qseaSet, "notAList"), "named list")
    expect_error(setParameters(qseaSet, list(1)), "named list")
})


test_that("setCNV() replaces the CNV ranges", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)
    skip_if_not(qsea::hasCNV(qseaSet))

    cnv <- qsea::getCNV(qseaSet)
    GenomicRanges::mcols(cnv)[] <- 0

    expect_equal(
        sum(as.matrix(GenomicRanges::mcols(qsea::getCNV(
            setCNV(qseaSet, cnv)
        )))),
        0
    )
})


test_that("the accessors round-trip through subsetQset()", {

    qseaSet <- cachedExampleQset(repl = 3, expSamplingDepth = 1000)
    samples <- qsea::getSampleNames(qseaSet)[1:2]

    newSet <- subsetQset(qseaSet, samplesToKeep = samples)

    expect_equal(qsea::getSampleNames(newSet), samples)
    expect_equal(colnames(qsea::getCounts(newSet)), samples)
    expect_equal(rownames(getLibraries(newSet, "file_name")), samples)
    expect_equal(rownames(qsea::getZygosity(newSet)), samples)
    expect_equal(
        rownames(getEnrichment(newSet, "parameters")), samples
    )
    expect_equal(colnames(getEnrichment(newSet, "factors")), samples)
})
