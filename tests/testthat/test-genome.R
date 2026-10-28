test_that("getMesaGenome() reads the mesa_genome option", {

    oldGenome <- getOption("mesa_genome")
    on.exit(options(mesa_genome = oldGenome), add = TRUE)

    options(mesa_genome = NULL)
    expect_equal(getMesaGenome(), "hg38")

    setMesaGenome("mm10")
    expect_equal(getMesaGenome(), "mm10")
})

test_that("getMesaAnnoDb() maps genomes to org packages", {

    oldGenome <- getOption("mesa_genome")
    on.exit(options(mesa_genome = oldGenome), add = TRUE)

    expect_equal(getMesaAnnoDb("hg38"), "org.Hs.eg.db")
    expect_equal(getMesaAnnoDb("hg19"), "org.Hs.eg.db")
    expect_equal(getMesaAnnoDb("mm10"), "org.Mm.eg.db")
    expect_error(getMesaAnnoDb("hg5"), "Unsupported genome: hg5")

    # NULL falls back to the mesa_genome option
    setMesaGenome("mm10")
    expect_equal(getMesaAnnoDb(), "org.Mm.eg.db")
})

test_that("getMesaTxDb() returns the TxDb for a genome", {

    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    oldGenome <- getOption("mesa_genome")
    on.exit(options(mesa_genome = oldGenome), add = TRUE)

    txdbGenome <- function(txdb) unique(GenomeInfoDb::genome(txdb))

    expect_s4_class(getMesaTxDb("hg38"), "TxDb")
    expect_equal(txdbGenome(getMesaTxDb("hg38")), "hg38")
    expect_equal(txdbGenome(getMesaTxDb("hg19")), "hg19")
    expect_equal(txdbGenome(getMesaTxDb("mm10")), "mm10")
    expect_error(getMesaTxDb("hg5"), "Unsupported genome: hg5")

    setMesaGenome("mm10")
    expect_equal(txdbGenome(getMesaTxDb()), "mm10")
})
