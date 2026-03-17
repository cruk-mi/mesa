#' Get current mesa genome setting
#' @return Character string of current genome setting (e.g., "hg38", "hg19", "mm10")
#' @export 
getMesaGenome <- function() {
  getOption("mesa_genome", "hg38")
}

#' Get TxDb for current or specified genome
#' @param genome Genome build, defaults to current setting
#' @return A TxDb object for the specified genome
#' @export
getMesaTxDb <- function(genome = NULL) {
  if (is.null(genome)) genome <- getMesaGenome()
  
  switch(genome,
         "hg38" = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
         "hg19" = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
         "mm10" = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
         stop("Unsupported genome: ", genome)
  )
}

#' Get annotation DB for current or specified genome
#' @param genome Genome build, defaults to current setting
#' @return Character string of annotation database name (e.g., "org.Hs.eg.db", "org.Mm.eg.db")
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