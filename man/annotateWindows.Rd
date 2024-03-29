% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qseaExtra.R
\name{annotateWindows}
\alias{annotateWindows}
\title{Add annotation onto a data table}
\usage{
annotateWindows(
  dataTable,
  genome = .getMesaGenome(),
  TxDb = .getMesaTxDb(),
  annoDb = .getMesaAnnoDb(),
  CpGislandsGR = NULL,
  FantomRegionsGR = NULL
)
}
\arguments{
\item{dataTable}{A data frame which can be coerced into a GRanges object or a GRanges object directly.}

\item{genome}{A genome string to set the rest of the parameters (currently only hg38/GRCh38 supported)}

\item{TxDb}{A TxDb database object (unquoted) to pass to ChIPseeker::annotatePeak}

\item{annoDb}{A string giving a Bioconductor annotation package, such as "org.Hs.eg.db"}

\item{CpGislandsGR}{A GRanges object giving locations of CpG islands}

\item{FantomRegionsGR}{A GRanges object giving Fantom enhancer regions}
}
\value{
A tibble with the data, augmented with ChIPseeker region location and CpG island information.
}
\description{
This function uses the ChIPseeker::annotatePeak function to determine the closest region to each genomic window provided.
Defaults to being hg38, unless
}
