% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyseDMRs.R
\name{summariseDMRsByContrast}
\alias{summariseDMRsByContrast}
\title{This function summarises the number of windows called as up/down in each comparison in a .
Internally transforms the data to a long format. Contrasts with no significant windows will be dropped!}
\usage{
summariseDMRsByContrast(
  DMRtable,
  FDRthres = 0.05,
  log2FCthres = 0,
  deltaBetaThres = 0
)
}
\arguments{
\item{DMRtable}{A data frame with multiple comparisons inside}

\item{FDRthres}{FDR threshold to apply to each comparison}

\item{log2FCthres}{A log2FC threshold to apply to each comparison (absolute)}

\item{deltaBetaThres}{A deltaBeta (change in average beta values) threshold to apply to each comparison (absolute)}
}
\description{
This function summarises the number of windows called as up/down in each comparison in a .
Internally transforms the data to a long format. Contrasts with no significant windows will be dropped!
}
