% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combineQsets.R
\name{combineQsetsList}
\alias{combineQsetsList}
\title{Combine two qseaSets}
\usage{
combineQsetsList(
  qseaSets,
  firstQset = NULL,
  dropDuplicates = TRUE,
  checkParams = TRUE,
  regionsToKeep = NULL
)
}
\arguments{
\item{qseaSets}{A list of qseaSet objects or strings to rds files with them in}

\item{firstQset}{An first qseaSet object to combine the rest into (optional)}

\item{dropDuplicates}{Whether to drop samples with the same name. Else renames them by adding "_Dup" to the end of the name.}

\item{checkParams}{Whether to verify that all the parameters are identical in both objects}

\item{regionsToKeep}{A GRanges object (or table coercible to one) to use to subset the samples. Helps keep RAM use down.}
}
\value{
A qseaSet object, containing all the samples from both qseaSet objects.
}
\description{
This function takes a list of qseaSets (or strings to read rds files from) and merges them together
}
