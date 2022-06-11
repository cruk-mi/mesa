% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotGeneHeatmap}
\alias{plotGeneHeatmap}
\title{This function takes a qseaSet and a gene, and plots the expression across the gene as a heatmap}
\usage{
plotGeneHeatmap(
  qseaSet,
  gene,
  normMethod = "beta",
  annotationCol = NULL,
  minEnrichment = 3,
  maxScale = 1,
  clusterNum = 2,
  annotationColors = NA,
  upstreamDist = 3000,
  scaleRows = FALSE,
  clusterCols = TRUE,
  downstreamDist = 1000
)
}
\arguments{
\item{qseaSet}{The qseaSet object.}

\item{gene}{A gene to plot. Currently as a gene id only, in a dubious way.}

\item{normMethod}{Whether to plot nrpm values or beta values.}

\item{annotationCol}{A data frame with annotations for the samples.}

\item{minEnrichment}{Minimum enrichment factor for beta values, will give NAs below this.}

\item{maxScale}{The maximum of the scale, not used when plotting beta values.}

\item{clusterNum}{A number of clusters to break the column dendrogram into.}

\item{annotationColors}{A list with the colours to use for the column legend, to pass to pheatmap}

\item{upstreamDist}{Number of basepairs upstream of the gene to include.}

\item{scaleRows}{Whether to scale the rows of the heatmap}

\item{downstreamDist}{Number of basepairs downstream of the gene to include.}
}
\value{
A qseaSet object with the sampleTable enhanced with the information on number of reads etc
}
\description{
This function takes a qseaSet and a gene, and plots the expression across the gene as a heatmap
}