#' Get current mesa genome setting
#'
#' Returns the genome identifier stored in the current R session for mesa's
#' annotation helpers.
#'
#' @return Character string of current genome setting (e.g., "hg38", "hg19",
#' "mm10")
#'
#' @examples
#' oldGenome <- getOption("mesa_genome")
#' options(mesa_genome = "hg19")
#' getMesaGenome()
#'
#' # Restore the previous genome setting
#' options(mesa_genome = oldGenome)
#' @export
getMesaGenome <- function() {
    getOption("mesa_genome", "hg38")
}

#' Get TxDb for current or specified genome
#'
#' @param genome Genome build, defaults to current setting
#' @return A TxDb object for the specified genome
#'
#' @examples
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'     getMesaTxDb("hg38")
#' }
#' @export
getMesaTxDb <- function(genome = NULL) {
    if (is.null(genome)) genome <- getMesaGenome()

    switch(genome,
        "hg38" = getExportedValue(
            "TxDb.Hsapiens.UCSC.hg38.knownGene",
            "TxDb.Hsapiens.UCSC.hg38.knownGene"
        ),
        "hg19" = getExportedValue(
            "TxDb.Hsapiens.UCSC.hg19.knownGene",
            "TxDb.Hsapiens.UCSC.hg19.knownGene"
        ),
        "mm10" = getExportedValue(
            "TxDb.Mmusculus.UCSC.mm10.knownGene",
            "TxDb.Mmusculus.UCSC.mm10.knownGene"
        ),
        stop("Unsupported genome: ", genome)
    )
}

#' Get annotation DB for current or specified genome
#'
#' @param genome Genome build, defaults to current setting
#' @return Character string of annotation database name (e.g., "org.Hs.eg.db",
#' "org.Mm.eg.db")
#'
#' @examples
#' getMesaAnnoDb("hg38")
#' getMesaAnnoDb("mm10")
#' @export
getMesaAnnoDb <- function(genome = NULL) {
    if (is.null(genome)) genome <- getMesaGenome()

    switch(genome,
        "hg38" = "org.Hs.eg.db",
        "hg19" = "org.Hs.eg.db",
        "mm10" = "org.Mm.eg.db",
        stop("Unsupported genome: ", genome)
    )
}
