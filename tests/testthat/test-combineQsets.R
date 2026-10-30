# Splits of exampleTumourNormal (10 samples, 819 windows) recombined in
# different ways.

allSamples <- exampleTumourNormal

tumours <- function() {
    filter(allSamples, tumour == "Tumour")
}
normals <- function() {
    filter(allSamples, tumour == "Normal")
}

tumourSamples <- qsea::getSampleNames(tumours())

test_that("combineQsets() renames clashing samples with _Dup", {

    expect_message(
        combined <- combineQsets(tumours(), tumours()),
        "adding _Dup to their name"
    )
    expect_equal(
        qsea::getSampleNames(combined),
        c(tumourSamples, paste0(tumourSamples, "_Dup"))
    )
    expect_equal(ncol(combined@count_matrix), 10)
    expect_equal(nrow(combined@zygosity), 10)
    expect_equal(nrow(combined@libraries$file_name), 10)
})

test_that("combineQsets(dropDuplicates = TRUE) drops shared samples", {

    # everything in the second set is already in the first
    expect_message(
        same <- combineQsets(
            exampleTumourNormal, exampleTumourNormal, dropDuplicates = TRUE
        ),
        "All names in common, returning first qseaSet"
    )
    expect_equal(qsea::getSampleNames(same),
        qsea::getSampleNames(exampleTumourNormal))

    # everything in the first set is already in the second
    expect_message(
        second <- combineQsets(
            tumours(), exampleTumourNormal, dropDuplicates = TRUE
        ),
        "All names in common, returning second qseaSet"
    )
    expect_equal(qsea::getSampleNames(second),
        qsea::getSampleNames(exampleTumourNormal))

    # partial overlap: the shared sample is kept once, from the first set
    overlap <- subsetQset(
        exampleTumourNormal, samplesToKeep = c("Colon1_T", "Colon1_N")
    )
    expect_message(
        partial <- combineQsets(tumours(), overlap, dropDuplicates = TRUE),
        "Dropping 1 common names"
    )
    expect_equal(qsea::getSampleNames(partial), c(tumourSamples, "Colon1_N"))
})

test_that("combineQsets() takes the intersection of differing regions", {

    keep <- qsea::getRegions(exampleTumourNormal)[1:400]
    reduced <- filterByOverlaps(normals(), keep)

    expect_message(
        combined <- combineQsets(tumours(), reduced),
        "Taking intersection of 400 regions"
    )
    expect_length(qsea::getRegions(combined), 400)
    expect_equal(nrow(combined@count_matrix), 400)
    expect_equal(ncol(combined@count_matrix), 10)
})

test_that("combineQsets(regionsToKeep) filters both sets first", {

    keep <- qsea::getRegions(exampleTumourNormal)[1:100]
    combined <- combineQsets(tumours(), normals(), regionsToKeep = keep)

    expect_length(qsea::getRegions(combined), 100)
    expect_equal(ncol(combined@count_matrix), 10)
})

test_that("combineQsets() warns when CpG densities differ", {

    other <- normals()
    other@regions$CpG_density <- other@regions$CpG_density + 1

    expect_warning(
        combined <- combineQsets(tumours(), other),
        "CpG densities differ between the two objects"
    )
    # the densities of the first set are kept
    expect_equal(
        qsea::getRegions(combined)$CpG_density,
        qsea::getRegions(tumours())$CpG_density
    )
})

test_that("combineQsets(checkParams = TRUE) stops on differing settings", {

    checked <- combineQsets(tumours(), normals(), checkParams = TRUE)
    expect_setequal(
        qsea::getSampleNames(checked), qsea::getSampleNames(allSamples)
    )

    other <- normals()
    other@parameters$minMapQual <- other@parameters$minMapQual + 10
    expect_error(
        combineQsets(tumours(), other, checkParams = TRUE),
        "Parameters are not the same!"
    )
    # the same difference is ignored by default
    ignored <- combineQsets(tumours(), other)
    expect_setequal(
        qsea::getSampleNames(ignored), qsea::getSampleNames(allSamples)
    )

    other <- normals()
    other@enrichment$extra <- 1
    expect_error(
        combineQsets(tumours(), other, checkParams = TRUE),
        "Enrichment entries are not the same"
    )

    other <- normals()
    other@enrichment$pattern_name <- "other"
    expect_error(
        combineQsets(tumours(), other, checkParams = TRUE),
        "Pattern names are different"
    )

    other <- normals()
    other@enrichment$density <- other@enrichment$density + 1
    expect_error(
        combineQsets(tumours(), other, checkParams = TRUE),
        "Enrichment density is different."
    )

    other <- normals()
    other@enrichment$n <- other@enrichment$n + 1
    expect_error(
        combineQsets(tumours(), other, checkParams = TRUE),
        "Enrichment n is different"
    )
})

test_that("combineQsets() needs qseaSets in both positions", {

    expect_error(
        combineQsets(1, tumours()),
        "Please supply a qseaSet object in the first position."
    )
    expect_error(
        combineQsets(tumours(), "notAnRds.txt"),
        "Please supply a qseaSet object in the second position."
    )
})

test_that("combineQsets() loads qseaSets from .rds paths", {

    first <- tempfile(fileext = ".rds")
    second <- tempfile(fileext = ".rds")
    on.exit(unlink(c(first, second)), add = TRUE)
    saveRDS(tumours(), first)
    saveRDS(normals(), second)

    expect_message(
        combined <- combineQsets(first, second),
        "Character string given, loading"
    )
    expect_setequal(
        qsea::getSampleNames(combined),
        qsea::getSampleNames(exampleTumourNormal)
    )
})

test_that("combineQsetsList() takes a first qseaSet or starts from the list", {

    sets <- list(a = tumours(), b = normals())

    expect_message(
        fromList <- combineQsetsList(sets),
        "using first element as initial qseaSet"
    )
    expect_setequal(
        qsea::getSampleNames(fromList),
        qsea::getSampleNames(exampleTumourNormal)
    )

    firstPath <- tempfile(fileext = ".rds")
    on.exit(unlink(firstPath), add = TRUE)
    saveRDS(tumours(), firstPath)
    expect_message(
        fromPath <- combineQsetsList(
            list(b = normals()), firstQset = firstPath
        ),
        "Character string given as firstQset"
    )
    expect_setequal(
        qsea::getSampleNames(fromPath),
        qsea::getSampleNames(exampleTumourNormal)
    )

    # a single set and no firstQset leaves nothing to start from
    expect_error(
        combineQsetsList(list(a = tumours())),
        "Please supply a qseaSet object in the first position."
    )
})
