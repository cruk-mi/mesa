% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilityFunctions.R
\name{setMesaTxDb}
\alias{setMesaTxDb}
\title{This function sets a default TxDb for the annotateWindows function to use.}
\usage{
setMesaTxDb(TxDb)
}
\arguments{
\item{TxDb}{A TxDb package (as a string, or the object itself)}
}
\value{
None
}
\description{
For instance "TxDb.Hsapiens.UCSC.hg38.knownGene" for human, or "TxDb.Mmusculus.UCSC.mm10.knownGene" for mouse.
}
\seealso{
The ChIPseeker function  used by annotateWindows is annotatePeak, \link[ChIPseeker]{annotatePeak}
}
