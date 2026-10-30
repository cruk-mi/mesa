# exampleTumourNormal has no low-CpG windows, so addNormalisation() cannot
# fit: every mix here uses renormalise = FALSE.

mixNReads <- 1e6

test_that("mixSamples() builds one mixed sample with default names", {

    set.seed(1)
    mixed <- suppressMessages(mixSamples(
        exampleTumourNormal, "Colon1_T", "Colon1_N", mixNReads, 0.3,
        onlyNew = TRUE, renormalise = FALSE
    ))
    newName <- "Mix_Colon1_T_Colon1_N_0.3"
    sampleTable <- qsea::getSampleTable(mixed)

    expect_equal(qsea::getSampleNames(mixed), newName)
    expect_equal(sampleTable$group, newName)
    expect_equal(sampleTable$sample1, "Colon1_T")
    expect_equal(sampleTable$sample2, "Colon1_N")
    expect_equal(sampleTable$prop1, 0.3)
    expect_equal(
        unname(mixed@libraries$file_name[newName, "valid_fragments"]),
        mixNReads
    )
    expect_equal(nrow(mixed@count_matrix), 819)
    expect_equal(colnames(mixed@count_matrix), newName)
})

test_that("mixSamples() draws reads from each sample in proportion", {

    set.seed(2)
    mixed <- suppressMessages(mixSamples(
        exampleTumourNormal, "Colon1_T", "Colon1_N", mixNReads, 0.3,
        onlyNew = TRUE, renormalise = FALSE
    ))
    libs <- exampleTumourNormal@libraries$file_name
    counts <- exampleTumourNormal@count_matrix
    frac <- colSums(counts[, c("Colon1_T", "Colon1_N")]) /
        libs[c("Colon1_T", "Colon1_N"), "valid_fragments"]
    expectedTotal <- floor(ceiling(mixNReads * 0.3) * frac[1]) +
        floor((mixNReads - ceiling(mixNReads * 0.3)) * frac[2])

    expect_equal(sum(mixed@count_matrix), unname(expectedTotal))
    # a mix cannot have more reads in a window than its two parents together
    expect_true(all(
        mixed@count_matrix[, 1] <=
            counts[, "Colon1_T"] + counts[, "Colon1_N"]
    ))
})

test_that("mixSamples() keeps the originals and honours newName/groupName", {

    set.seed(3)
    mixed <- suppressMessages(mixSamples(
        exampleTumourNormal, "Lung1_T", "Lung1_N", mixNReads, 0.5,
        newName = "mixA", groupName = "mixGroup", renormalise = FALSE
    ))
    sampleTable <- qsea::getSampleTable(mixed)

    expect_equal(
        qsea::getSampleNames(mixed),
        c(qsea::getSampleNames(exampleTumourNormal), "mixA")
    )
    expect_equal(sampleTable[sampleTable$sample_name == "mixA", "group"],
        "mixGroup")
    expect_equal(sampleTable[sampleTable$sample_name == "Lung1_T", "group"],
        "Lung1_T")
    # columns that only the originals have are NA for the mix
    expect_true(is.na(sampleTable["mixA", "patient"]))
    expect_equal(sampleTable["mixA", "prop1"], 0.5)
})

test_that("mixSamples() mixes the CNV and zygosity by proportion", {

    set.seed(4)
    prop <- 0.25
    mixed <- suppressMessages(mixSamples(
        exampleTumourNormal, "Lung1_T", "Lung1_N", mixNReads, prop,
        newName = "mixCNV", onlyNew = TRUE, renormalise = FALSE
    ))
    cnv <- as.matrix(GenomicRanges::mcols(qsea::getCNV(exampleTumourNormal)))
    zyg <- exampleTumourNormal@zygosity

    expect_equal(
        unname(as.matrix(GenomicRanges::mcols(qsea::getCNV(mixed)))[, 1]),
        unname(prop * cnv[, "Lung1_T"] + (1 - prop) * cnv[, "Lung1_N"])
    )
    expect_equal(
        unname(mixed@zygosity["mixCNV", ]),
        unname(prop * zyg["Lung1_T", ] + (1 - prop) * zyg["Lung1_N", ])
    )
})

test_that("mixSamples() stops on bad input", {

    expect_error(
        mixSamples(exampleTumourNormal, "Colon1_T", "nope", mixNReads, 0.3),
        "Sample(s) not found in qseaSet: nope", fixed = TRUE
    )
    expect_error(
        suppressMessages(mixSamples(
            exampleTumourNormal, "Colon1_T", "Colon1_N", mixNReads, 1.2,
            renormalise = FALSE
        )),
        "Proportion must be a number between 0 and 1, not 1.2."
    )
    expect_error(
        suppressMessages(mixSamples(
            exampleTumourNormal, "Colon1_T", "Colon1_N", 1e9, 0.5,
            renormalise = FALSE
        )),
        "Not enough reads in Colon1_T"
    )
})

test_that("mixThreeQsetSamples() stops on bad input", {

    # Only the error paths are tested: its success path currently fails
    # (parked in the PR).
    mix3 <- function(...) {
        suppressMessages(mixThreeQsetSamples(
            exampleTumourNormal, "Colon1_T", "Colon1_N", "Lung1_T",
            nReadsTotal = mixNReads, renormalise = FALSE, ...
        ))
    }

    expect_error(
        mixThreeQsetSamples(
            exampleTumourNormal, "Colon1_T", "Colon1_N", "nope",
            mixNReads, 0.2, 0.3
        ),
        "Sample(s) not found in qseaSet: nope", fixed = TRUE
    )
    expect_error(
        mix3(proportion1 = 1.5, proportion2 = 0.1),
        "Proportion1 must be a number between 0 and 1, not 1.5."
    )
    expect_error(
        mix3(proportion1 = 0.2, proportion2 = -0.1),
        "Proportion2 must be a number between 0 and 1, not -0.1."
    )
    expect_error(
        mix3(proportion1 = 0.6, proportion2 = 0.6),
        "Sum of proportions must be less than 1, not 1.2."
    )
})
