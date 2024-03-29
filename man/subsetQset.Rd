% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alterQset.R
\name{subsetQset}
\alias{subsetQset}
\title{Subset a qseaSet}
\usage{
subsetQset(qseaSet, samplesToKeep = NULL, samplesToDrop = NULL)
}
\arguments{
\item{qseaSet}{The original qseaSet object.}

\item{samplesToKeep}{Which samples to keep. Only one of samplesToKeep or samplesToDrop should be specified.}

\item{samplesToDrop}{Which samples to drop.}
}
\value{
A qseaSet object, with only the selected samples inside it.
}
\description{
This function takes a qseaSet and keeps only a subset of the samples from it.
}
\examples{
subsetQset(exampleTumourNormal, samplesToKeep = c('Colon1_T','Colon1_N'))
}
